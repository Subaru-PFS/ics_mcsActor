/* Spots of a frame: level and noise of its columns, detection, windowed centroids and fluxes. */
#include <math.h>
#include <pthread.h>
#include <stdlib.h>
#include "spotFinder.h"

/* Bins, 1 ADU wide, of the histogram a column is reduced with. */
#define NBINS 512

/* Fraction of a Gaussian below one standard deviation under its mean. */
#define LOWER_FRACTION 0.15865525393145707
/* Side, in pixels, of the blocks whose medians make the smooth background of a frame. */
#define MESH_SIZE 256
/* Step, in rows, of the pixels of a block whose median is taken: every column, every 4th row, about 16000
   pixels, the median good to 0.1 ADU at the noise of a frame. Every column is kept, so that each block
   holds the same mixture of column levels. */
#define MESH_STEP 4

/* Lowest detection threshold, ADU: a frame without noise is not searched at zero. */
#define MIN_THRESHOLD 1.0

/* Pixels above the threshold beyond which a frame is not searched; a frame of the camera has about 150 000. */
#define MAX_BRIGHT_PIXELS 3000000

/* Spots beyond which a frame is given up. */
#define MAX_SPOTS 10000

#define MAX_THREADS 64

/* Change of position, in pixels, under which an axis of the windowed centroid has converged. */
#define PRECISION 1e-4

/*------------------------------------------------------------------------------------------*/

typedef void (*rangeFunction)(void *context, int first, int last);

typedef struct rangeJob
{
  rangeFunction function;
  void *context;
  int first;
  int last;
} rangeJob;

static void *runRange(void *arg)
{
  rangeJob *job = (rangeJob *)arg;
  job->function(job->context, job->first, job->last);
  return NULL;
}

static void overRanges(rangeFunction function, void *context, int n, int nThreads)
{
  /* call function(context, first, last) on consecutive pieces of 0..n, each in its own thread.
     A piece whose thread cannot be created runs in the calling one. */

  pthread_t threads[MAX_THREADS];
  rangeJob jobs[MAX_THREADS];
  int started[MAX_THREADS];
  int k;

  if(n < 1)
    {
      return;
    }
  if(nThreads > MAX_THREADS)
    {
      nThreads = MAX_THREADS;
    }
  if(nThreads > n)
    {
      nThreads = n;
    }
  if(nThreads < 1)
    {
      nThreads = 1;
    }

  for(k = 0; k < nThreads; k++)
    {
      jobs[k].function = function;
      jobs[k].context = context;
      jobs[k].first = (int)((long)n * k / nThreads);
      jobs[k].last = (int)((long)n * (k + 1) / nThreads);
    }
  for(k = 1; k < nThreads; k++)
    {
      started[k] = (pthread_create(&threads[k], NULL, runRange, &jobs[k]) == 0);
      if(!started[k])
	{
	  runRange(&jobs[k]);
	}
    }
  runRange(&jobs[0]);
  for(k = 1; k < nThreads; k++)
    {
      if(started[k])
	{
	  pthread_join(threads[k], NULL);
	}
    }
}

/*------------------------------------------------------------------------------------------*/

static int compareInt(const void *a, const void *b)
{
  int x = *(const int *)a, y = *(const int *)b;
  return (x > y) - (x < y);
}

static int compareDouble(const void *a, const void *b)
{
  double x = *(const double *)a, y = *(const double *)b;
  return (x > y) - (x < y);
}

static int sampleMedian(const int *image, long n)
{
  /* median of a regular sample of about 4000 pixels */
  enum { NSAMPLE = 4001 };
  int sample[NSAMPLE];
  long step = n / NSAMPLE;
  int count = 0;
  long i;

  if(step < 1)
    {
      step = 1;
    }
  for(i = step / 2; i < n && count < NSAMPLE; i += step)
    {
      sample[count++] = image[i];
    }
  qsort(sample, count, sizeof(int), compareInt);
  return sample[count / 2];
}

static double median(const double *values, int n)
{
  /* as numpy.median; values are left as they are. NaN when memory is short. */
  double *sorted = malloc(sizeof(double) * n);
  double result;
  int k;

  if(sorted == NULL)
    {
      return NAN;
    }
  for(k = 0; k < n; k++)
    {
      sorted[k] = values[k];
    }
  qsort(sorted, n, sizeof(double), compareDouble);
  result = 0.5 * (sorted[(n - 1) / 2] + sorted[n / 2]);
  free(sorted);
  return result;
}

static int sortedColumn(const int *image, int n_x, int n_y, int x, const int *ranks, double *values)
{
  /* the four order statistics of column x, by sorting it; returns 0 when memory is short */
  int *column = malloc(sizeof(int) * n_y);
  int y, k;

  if(column == NULL)
    {
      return 0;
    }
  for(y = 0; y < n_y; y++)
    {
      column[y] = image[(long)y * n_x + x];
    }
  qsort(column, n_y, sizeof(int), compareInt);
  for(k = 0; k < 4; k++)
    {
      values[k] = column[ranks[k]];
    }
  free(column);
  return 1;
}

static int histogramColumn(const unsigned short *histogram, int below, int low, const int *ranks, double *values)
{
  /* the four order statistics of a column from its histogram; ranks ascend.
     Returns 0 when one of them lies outside the histogram. */
  int cumulative = below, bin = 0, k;

  for(k = 0; k < 4; k++)
    {
      if(ranks[k] < below)
	{
	  return 0;
	}
      while(bin < NBINS && cumulative + histogram[bin] <= ranks[k])
	{
	  cumulative += histogram[bin];
	  bin++;
	}
      if(bin == NBINS)
	{
	  return 0;
	}
      values[k] = low + bin;
    }
  return 1;
}

typedef struct columnContext
{
  const int *image;
  int n_x;
  int n_y;
  int low;            /* first value of the histograms */
  int ranks[4];       /* order statistics around the lower percentile, then around the median */
  double fraction;    /* weight of the second one in the lower percentile */
  double *levels;     /* per column: median */
  double *lower;      /* per column: median minus lower percentile */
  unsigned char *sorted;  /* per column: 1 when reduced by sorting, 2 when it could not be reduced */
} columnContext;

static void reduceColumns(void *context, int first, int last)
{
  /* levels and lower widths of columns first..last-1, the image being read in row order */

  const columnContext *c = (const columnContext *)context;
  const int width = last - first;
  unsigned short *histogram = calloc((size_t)width * NBINS, sizeof(unsigned short));
  unsigned short *below = calloc(width, sizeof(unsigned short));
  double values[4];
  int x, y;

  if(histogram != NULL && below != NULL)
    {
      for(y = 0; y < c->n_y; y++)
	{
	  const int *row = c->image + (long)y * c->n_x + first;
	  for(x = 0; x < width; x++)
	    {
	      /* widened, so that a value far from the level cannot wrap into the window */
	      const long offset = (long)row[x] - c->low;
	      if(offset < 0)
		{
		  below[x]++;
		}
	      else if(offset < NBINS)
		{
		  histogram[(size_t)x * NBINS + offset]++;
		}
	    }
	}
    }

  for(x = 0; x < width; x++)
    {
      /* the ranks ascend, as histogramColumn needs, from three rows on */
      if(histogram == NULL || below == NULL || c->n_y < 3
	 || !histogramColumn(histogram + (size_t)x * NBINS, below[x], c->low, c->ranks, values))
	{
	  if(!sortedColumn(c->image, c->n_x, c->n_y, first + x, c->ranks, values))
	    {
	      c->sorted[first + x] = 2;
	      continue;
	    }
	  c->sorted[first + x] = 1;
	}
      c->levels[first + x] = 0.5 * (values[2] + values[3]);
      c->lower[first + x] = c->levels[first + x] - (values[0] + c->fraction * (values[1] - values[0]));
    }

  free(histogram);
  free(below);
}

int columnLevels(const int *image, int n_x, int n_y, int nThreads, double *levels, double *noise)
{
  /*
    Median of every column, and the noise of the image.

    input
      image: image, n_y rows of n_x pixels
      nThreads: threads used at most
    output
      levels[n_x]: median over the rows of each column, as numpy.median(image, axis=0)
      noise: root mean square, over the columns, of the median minus the 15.87th percentile (linearly
             interpolated between order statistics). Light can only raise a pixel, so the lower side of a
	     column's distribution holds its noise alone.

    Returns the number of columns sorted because an order statistic fell outside their histogram,
    -1 when the image is empty or has more rows than a histogram bin can count, -2 when memory is short.

    The histograms cover NBINS ADU around the level of the image.
  */

  columnContext context;
  const double position = LOWER_FRACTION * (n_y - 1);
  double sum = 0.;
  int nSorted = 0;
  int x;

  if(n_x < 1 || n_y < 1 || n_y > 65535)
    {
      return -1;
    }

  context.image = image;
  context.n_x = n_x;
  context.n_y = n_y;
  context.low = sampleMedian(image, (long)n_x * n_y) - NBINS / 2;
  context.ranks[0] = (int)floor(position);
  context.ranks[1] = (int)ceil(position);
  context.ranks[2] = (n_y - 1) / 2;
  context.ranks[3] = n_y / 2;
  context.fraction = position - floor(position);
  context.levels = levels;
  context.lower = malloc(sizeof(double) * n_x);
  context.sorted = calloc(n_x, sizeof(unsigned char));
  if(context.lower == NULL || context.sorted == NULL)
    {
      free(context.lower);
      free(context.sorted);
      return -2;
    }

  overRanges(reduceColumns, &context, n_x, nThreads);

  for(x = 0; x < n_x; x++)
    {
      if(context.sorted[x] == 2)
	{
	  nSorted = -2;
	  break;
	}
      sum += context.lower[x] * context.lower[x];
      nSorted += context.sorted[x];
    }
  *noise = sqrt(sum / n_x);

  free(context.lower);
  free(context.sorted);
  return nSorted;
}

/*------------------------------------------------------------------------------------------*/

static int selectKth(int *values, long n, long k)
{
  /* the k-th smallest of values (0 based), which are reordered */
  long low = 0, high = n - 1;

  while(high > low)
    {
      const int pivot = values[low + (high - low) / 2];
      long i = low, j = high;
      while(i <= j)
	{
	  while(values[i] < pivot)
	    {
	      i++;
	    }
	  while(values[j] > pivot)
	    {
	      j--;
	    }
	  if(i <= j)
	    {
	      const int t = values[i];
	      values[i++] = values[j];
	      values[j--] = t;
	    }
	}
      if(k <= j)
	{
	  high = j;
	}
      else if(k >= i)
	{
	  low = i;
	}
      else
	{
	  break;
	}
    }
  return values[k];
}

typedef struct meshContext
{
  int *image;
  int n_x;
  int n_y;
  int nBlocksX;
  int nBlocksY;
  double *mesh;           /* nBlocksY rows of nBlocksX block medians */
  double level;           /* median of the block medians */
  const int *columnBlock; /* per column: the block on the left of the column, for the interpolation */
  const double *columnWeight;  /* per column: weight of the block on its right */
  const int *rowBlock;
  const double *rowWeight;
} meshContext;

static int blockStart(int block, int nBlocks, int n)
{
  return (int)((long)block * n / nBlocks);
}

static void meshMedians(void *context, int first, int last)
{
  /* the medians of the blocks of block rows first..last-1, of every MESH_STEP-th row; NAN when memory is
     short */
  const meshContext *c = (const meshContext *)context;
  const int width = c->n_x / c->nBlocksX + 1, height = c->n_y / c->nBlocksY / MESH_STEP + 2;
  int *values = malloc(sizeof(int) * (size_t)width * height);
  int i, j, x, y;

  for(j = first; j < last; j++)
    {
      const int y0 = blockStart(j, c->nBlocksY, c->n_y), y1 = blockStart(j + 1, c->nBlocksY, c->n_y);
      for(i = 0; i < c->nBlocksX; i++)
	{
	  const int x0 = blockStart(i, c->nBlocksX, c->n_x), x1 = blockStart(i + 1, c->nBlocksX, c->n_x);
	  long n = 0;
	  int lowerMiddle, upperMiddle;
	  if(values == NULL)
	    {
	      c->mesh[j * c->nBlocksX + i] = NAN;
	      continue;
	    }
	  for(y = y0; y < y1; y += MESH_STEP)
	    {
	      const int *row = c->image + (long)y * c->n_x;
	      for(x = x0; x < x1; x++)
		{
		  values[n++] = row[x];
		}
	    }
	  /* as numpy.median: the mean of the two middle values of an even number */
	  upperMiddle = selectKth(values, n, n / 2);
	  lowerMiddle = upperMiddle;
	  if(n % 2 == 0)
	    {
	      long k;
	      lowerMiddle = values[0];
	      for(k = 1; k < n / 2; k++)
		{
		  if(values[k] > lowerMiddle)
		    {
		      lowerMiddle = values[k];
		    }
		}
	    }
	  c->mesh[j * c->nBlocksX + i] = 0.5 * ((double)lowerMiddle + upperMiddle);
	}
    }
  free(values);
}

static void interpolationAxis(int n, int nBlocks, int *block, double *weight)
{
  /* for each pixel of an axis of n pixels in nBlocks blocks: the block whose centre is at or before it
     (the first one before the first centre, the one before last after the last), and the weight of the next
     block, extrapolated linearly beyond the first and last centres */
  int pixel, b = 0;

  for(pixel = 0; pixel < n; pixel++)
    {
      double centre, next;
      while(b < nBlocks - 2 && 0.5 * (blockStart(b + 1, nBlocks, n) + blockStart(b + 2, nBlocks, n) - 1) <= pixel)
	{
	  b++;
	}
      if(nBlocks == 1)
	{
	  block[pixel] = 0;
	  weight[pixel] = 0.;
	  continue;
	}
      centre = 0.5 * (blockStart(b, nBlocks, n) + blockStart(b + 1, nBlocks, n) - 1);
      next = 0.5 * (blockStart(b + 1, nBlocks, n) + blockStart(b + 2, nBlocks, n) - 1);
      block[pixel] = b;
      weight[pixel] = (pixel - centre) / (next - centre);
    }
}

static double meshValue(const meshContext *c, int x, int y)
{
  /* the smooth background at a pixel, interpolated bilinearly between block centres */
  const int i = c->columnBlock[x], j = c->rowBlock[y];
  const int i1 = c->nBlocksX > 1 ? i + 1 : i, j1 = c->nBlocksY > 1 ? j + 1 : j;
  const double wx = c->columnWeight[x], wy = c->rowWeight[y];
  const double *m = c->mesh;
  return (1. - wy) * ((1. - wx) * m[j * c->nBlocksX + i] + wx * m[j * c->nBlocksX + i1])
    + wy * ((1. - wx) * m[j1 * c->nBlocksX + i] + wx * m[j1 * c->nBlocksX + i1]);
}

static int meshDeviation(const meshContext *c, int x, int y)
{
  /* the smooth background at a pixel less the median of the blocks, rounded: 0 wherever the background is
     flat to half an ADU, so that a flat frame is left as it is */
  return (int)floor(meshValue(c, x, y) - c->level + 0.5);
}

static void subtractMesh(void *context, int first, int last)
{
  /* subtract the deviation of the smooth background from rows first..last-1, as meshDeviation gives it,
     interpolating each row between the block rows once */
  const meshContext *c = (const meshContext *)context;
  double *along = malloc(sizeof(double) * c->nBlocksX);
  int i, x, y;

  for(y = first; y < last; y++)
    {
      int *row = c->image + (long)y * c->n_x;
      const int j = c->rowBlock[y], j1 = c->nBlocksY > 1 ? j + 1 : j;
      const double wy = c->rowWeight[y];
      if(along == NULL)
	{
	  for(x = 0; x < c->n_x; x++)
	    {
	      row[x] -= meshDeviation(c, x, y);
	    }
	  continue;
	}
      for(i = 0; i < c->nBlocksX; i++)
	{
	  along[i] = (1. - wy) * c->mesh[j * c->nBlocksX + i] + wy * c->mesh[j1 * c->nBlocksX + i];
	}
      for(x = 0; x < c->n_x; x++)
	{
	  const int b = c->columnBlock[x], b1 = c->nBlocksX > 1 ? b + 1 : b;
	  const double wx = c->columnWeight[x];
	  row[x] -= (int)floor((1. - wx) * along[b] + wx * along[b1] - c->level + 0.5);
	}
    }
  free(along);
}

/*------------------------------------------------------------------------------------------*/

double apertureSum(const int *image, int n_x, int n_y, double x, double y, double radius, int *nPixels)
{
  /*
    Sum of the pixels whose centre lies within radius of (x, y), pixel (i, j) being centred on x = i, y = j,
    and their number. The image is taken to have no background. NaN, and no pixel, when the position is
    not on the image or the aperture leaves it.
  */

  const double radius2 = radius * radius;
  const int reach = (int)ceil(radius);
  double sum = 0.;
  int count = 0;
  int i, j;

  *nPixels = 0;
  if(!(x >= 0 && x <= n_x - 1 && y >= 0 && y <= n_y - 1))
    {
      return NAN;
    }

  const int i0 = (int)floor(x + 0.5);
  const int j0 = (int)floor(y + 0.5);

  if(i0 - reach < 0 || i0 + reach >= n_x || j0 - reach < 0 || j0 + reach >= n_y)
    {
      return NAN;
    }

  for(j = j0 - reach; j <= j0 + reach; j++)
    {
      const int *row = image + (long)j * n_x;
      const double dy2 = (j - y) * (j - y);
      for(i = i0 - reach; i <= i0 + reach; i++)
	{
	  if((i - x) * (i - x) + dy2 <= radius2)
	    {
	      sum += row[i];
	      count++;
	    }
	}
    }
  *nPixels = count;
  return sum;
}

/*------------------------------------------------------------------------------------------*/

typedef struct region
{
  double x;     /* centroid, first of the pixels above the threshold, then windowed */
  double y;
  double x2;
  double y2;
  double xy;
  double peak;
} region;

static int floodRegion(const int *image, unsigned char *mask, int n_x, int n_y, long start, int threshold,
		       long **stack, long *capacity, region *found)
{
  /*
    The pixels at or above the threshold connected, side to side, to the pixel start, which is above it.
    They are marked in mask. found gets their intensity-weighted centroid and second moments and their
    highest value. Returns their number, or -1 when the stack could not grow.
  */

  const int x0 = (int)(start % n_x), y0 = (int)(start / n_x);
  double total = 0., sumX = 0., sumY = 0., sumX2 = 0., sumY2 = 0., sumXY = 0.;
  int peak = image[start];
  long depth = 0;
  int nPixels = 0;
  int k;

  mask[start] = 1;
  (*stack)[depth++] = start;

  while(depth > 0)
    {
      const long pixel = (*stack)[--depth];
      const int i = (int)(pixel % n_x), j = (int)(pixel / n_x);
      const double value = image[pixel];
      const double dx = i - x0, dy = j - y0;
      const long neighbours[4] = {i + 1 < n_x ? pixel + 1 : -1, i > 0 ? pixel - 1 : -1,
				  j + 1 < n_y ? pixel + n_x : -1, j > 0 ? pixel - n_x : -1};

      total += value;
      sumX += value * dx;
      sumY += value * dy;
      sumX2 += value * dx * dx;
      sumY2 += value * dy * dy;
      sumXY += value * dx * dy;
      if(image[pixel] > peak)
	{
	  peak = image[pixel];
	}
      nPixels++;

      if(depth + 4 > *capacity)
	{
	  long *grown = realloc(*stack, sizeof(long) * 2 * (*capacity));
	  if(grown == NULL)
	    {
	      return -1;
	    }
	  *stack = grown;
	  *capacity *= 2;
	}
      for(k = 0; k < 4; k++)
	{
	  const long next = neighbours[k];
	  if(next >= 0 && !mask[next] && image[next] >= threshold)
	    {
	      mask[next] = 1;
	      (*stack)[depth++] = next;
	    }
	}
    }

  found->x = sumX / total + x0;
  found->y = sumY / total + y0;
  found->x2 = sumX2 / total - (sumX / total) * (sumX / total);
  found->y2 = sumY2 / total - (sumY / total) * (sumY / total);
  found->xy = sumXY / total - (sumX / total) * (sumY / total);
  found->peak = peak;
  return nPixels;
}

static double newtonGain(double window2, double variance)
{
  /*
    Factor from the weighted mean offset of the pixels to the offset of the spot, along one axis.

    For a Gaussian spot of variance s2 seen through a Gaussian window of variance window2, the weighted
    mean offset is the spot's offset times window2 / (window2 + s2), and the weighted variance is
    window2 s2 / (window2 + s2); the factor (window2 + s2) / window2 is then window2 / (window2 -
    variance). It is kept within 1 and 2, the constant factor of SExtractor, which is right only when the
    spot is as wide as the window. The point the centroid converges to does not depend on it.
  */

  double gain;

  if(!(variance < window2))
    {
      return 2.;
    }
  gain = window2 / (window2 - variance);
  if(gain < 1.)
    {
      return 1.;
    }
  return gain > 2. ? 2. : gain;
}

static void windowedCentroid(const int *image, int n_x, int n_y, const spotParams *params, double *x, double *y)
{
  /*
    Windowed centroid of the spot at (x, y), after SExtractor's windowed positions: the pixels within
    boxCent of the current position are weighted by a Gaussian of full width fwhmx for x and fwhmy for y,
    and the position moves by their weighted mean offset times newtonGain until both axes move by less
    than PRECISION, or maxIt iterations. The position is left where the weights no longer sum to a
    positive value.
  */

  const double sigmaX = params->fwhmx / sqrt(8. * log(2.));
  const double sigmaY = params->fwhmy / sqrt(8. * log(2.));
  const double box = params->boxCent;
  double xc = *x, yc = *y;
  double dx = 2 * PRECISION, dy = 2 * PRECISION;
  int iteration = 0;
  int i, j;

  while((dx > PRECISION || dy > PRECISION) && iteration < params->maxIt)
    {
      double sumX = 0., sumY = 0., sumX2 = 0., sumY2 = 0., weightX = 0., weightY = 0.;
      double meanX, meanY, stepX, stepY;
      int iMin = (int)floor(xc - box), iMax = (int)ceil(xc + box);
      int jMin = (int)floor(yc - box), jMax = (int)ceil(yc + box);

      if(iMin < 0)
	{
	  iMin = 0;
	}
      if(iMax > n_x - 1)
	{
	  iMax = n_x - 1;
	}
      if(jMin < 0)
	{
	  jMin = 0;
	}
      if(jMax > n_y - 1)
	{
	  jMax = n_y - 1;
	}

      {
	/* the Gaussians are separable: exp(-r^2 / 2s^2) = exp(-dx^2 / 2s^2) exp(-dy^2 / 2s^2) */
	const int nI = iMax - iMin + 1, nJ = jMax - jMin + 1;
	double offsetI[nI > 0 ? nI : 1], windowXofI[nI > 0 ? nI : 1], windowYofI[nI > 0 ? nI : 1];
	double offsetJ[nJ > 0 ? nJ : 1], windowXofJ[nJ > 0 ? nJ : 1], windowYofJ[nJ > 0 ? nJ : 1];

	for(i = iMin; i <= iMax; i++)
	  {
	    const double d = i - xc;
	    offsetI[i - iMin] = d;
	    windowXofI[i - iMin] = exp(-d * d / (2. * sigmaX * sigmaX));
	    windowYofI[i - iMin] = exp(-d * d / (2. * sigmaY * sigmaY));
	  }
	for(j = jMin; j <= jMax; j++)
	  {
	    const double d = j - yc;
	    offsetJ[j - jMin] = d;
	    windowXofJ[j - jMin] = exp(-d * d / (2. * sigmaX * sigmaX));
	    windowYofJ[j - jMin] = exp(-d * d / (2. * sigmaY * sigmaY));
	  }

	for(j = jMin; j <= jMax; j++)
	  {
	    const int *row = image + (long)j * n_x;
	    const double dy2 = offsetJ[j - jMin] * offsetJ[j - jMin];
	    for(i = iMin; i <= iMax; i++)
	      {
		const double di = offsetI[i - iMin];
		if(di * di + dy2 < box * box)
		  {
		    const double wx = windowXofI[i - iMin] * windowXofJ[j - jMin] * row[i];
		    const double wy = windowYofI[i - iMin] * windowYofJ[j - jMin] * row[i];
		    sumX += wx * di;
		    sumY += wy * offsetJ[j - jMin];
		    sumX2 += wx * di * di;
		    sumY2 += wy * dy2;
		    weightX += wx;
		    weightY += wy;
		  }
	      }
	  }
      }

      if(!(weightX > 0) || !(weightY > 0))
	{
	  break;
	}

      meanX = sumX / weightX;
      meanY = sumY / weightY;
      stepX = newtonGain(sigmaX * sigmaX, sumX2 / weightX - meanX * meanX) * meanX;
      stepY = newtonGain(sigmaY * sigmaY, sumY2 / weightY - meanY * meanY) * meanY;
      dx = fabs(stepX);
      dy = fabs(stepY);
      xc += stepX;
      yc += stepY;
      iteration++;

      if(!(xc >= 0 && xc <= n_x - 1 && yc >= 0 && yc <= n_y - 1))
	{
	  break;
	}
    }

  *x = xc;
  *y = yc;
}

typedef struct rankedY
{
  double y;
  int index;
} rankedY;

static int compareRankedY(const void *a, const void *b)
{
  /* by y, then by index, so that the order does not depend on qsort */
  const rankedY *p = (const rankedY *)a, *q = (const rankedY *)b;
  if(p->y != q->y)
    {
      return (p->y > q->y) - (p->y < q->y);
    }
  return (p->index > q->index) - (p->index < q->index);
}

static int columnOf(double x, int n_x)
{
  /* column of the pixel centred nearest to x, within the image */
  const double column = floor(x + 0.5);

  if(!(column >= 0))
    {
      return 0;
    }
  if(column > n_x - 1)
    {
      return n_x - 1;
    }
  return (int)column;
}

typedef struct centroidContext
{
  const int *image;
  int n_x;
  int n_y;
  const spotParams *params;
  region *regions;
} centroidContext;

static void centroidRegions(void *context, int first, int last)
{
  const centroidContext *c = (const centroidContext *)context;
  int k;

  for(k = first; k < last; k++)
    {
      windowedCentroid(c->image, c->n_x, c->n_y, c->params, &c->regions[k].x, &c->regions[k].y);
    }
}

spot *findSpots(int *image, int n_x, int n_y, const spotParams *params, frameLevels *levels, int *nSpots)
{
  /*
    The spots of a frame.

    The background of the frame is subtracted from the image, in place, in two steps: a smooth background,
    interpolated bilinearly between the medians of blocks of about MESH_SIZE pixels a side (of every
    MESH_STEP-th row), by its deviation
    from the median of the blocks rounded to an integer (see meshDeviation), then the level of each column,
    the median of what is left. The noise is that of what is left after the smooth background (see
    columnLevels): structure the background cannot follow, such as front-lit hardware, raises it, and the
    threshold with it.
    A spot is a set of
    pixels connected side to side, at or above the threshold and one of them above it, nmin pixels at
    least, whose centroid lies further than three times the larger of boxFind and boxCent from the edges
    of the image. The threshold is findSigma times the noise (see columnLevels), MIN_THRESHOLD at least,
    truncated to an integer. Of spots whose windowed centroids are closer than boxFind + 1 pixels, the one
    with the highest peak is kept.

    input
      image: frame, n_y rows of n_x pixels, of values within +-2^30; left with the background subtracted
      params: see spotParams
    output
      levels: background, noise and threshold of the frame
      nSpots: number of spots, or a negative status (see spotFinder.h)

    Returns the spots, allocated with malloc, in the order their first pixel comes in the image;
    NULL when nSpots is negative.
  */

  const long nPixels = (long)n_x * n_y;
  const double margin = 3. * (params->boxFind > params->boxCent ? params->boxFind : params->boxCent);
  const double closest = params->boxFind + 1;
  double *columns = NULL;
  int *rounded = NULL;
  long *bright = NULL, *stack = NULL;
  unsigned char *mask = NULL, *dropped = NULL;
  rankedY *byY = NULL;
  int *rankOf = NULL, *near = NULL;
  region *regions = NULL;
  spot *spots = NULL;
  centroidContext context;
  meshContext mesh;
  long nBright = 0, capacity = 4096, k;
  int nRegions = 0, nKept = 0;
  int status, threshold, x, y, a, b;

  *nSpots = SPOTS_NO_MEMORY;
  if(n_x < 1 || n_y < 1)
    {
      *nSpots = SPOTS_BAD_IMAGE;
      return NULL;
    }
  mesh.image = image;
  mesh.n_x = n_x;
  mesh.n_y = n_y;
  mesh.nBlocksX = (n_x + MESH_SIZE / 2) / MESH_SIZE > 0 ? (n_x + MESH_SIZE / 2) / MESH_SIZE : 1;
  mesh.nBlocksY = (n_y + MESH_SIZE / 2) / MESH_SIZE > 0 ? (n_y + MESH_SIZE / 2) / MESH_SIZE : 1;
  mesh.mesh = malloc(sizeof(double) * mesh.nBlocksX * mesh.nBlocksY);
  mesh.columnBlock = malloc(sizeof(int) * n_x);
  mesh.columnWeight = malloc(sizeof(double) * n_x);
  mesh.rowBlock = malloc(sizeof(int) * n_y);
  mesh.rowWeight = malloc(sizeof(double) * n_y);

  columns = malloc(sizeof(double) * n_x);
  rounded = malloc(sizeof(int) * n_x);
  bright = malloc(sizeof(long) * MAX_BRIGHT_PIXELS);
  mask = calloc(nPixels, sizeof(unsigned char));
  stack = malloc(sizeof(long) * capacity);
  regions = malloc(sizeof(region) * (MAX_SPOTS + 1));
  if(columns == NULL || rounded == NULL || bright == NULL || mask == NULL || stack == NULL || regions == NULL
     || mesh.mesh == NULL || mesh.columnBlock == NULL || mesh.columnWeight == NULL || mesh.rowBlock == NULL
     || mesh.rowWeight == NULL)
    {
      goto done;
    }

  /* the smooth background, then the level of each column and the noise in what is left */
  overRanges(meshMedians, &mesh, mesh.nBlocksY, params->nThreads);
  for(k = 0; k < (long)mesh.nBlocksX * mesh.nBlocksY; k++)
    {
      if(isnan(mesh.mesh[k]))
	{
	  goto done;
	}
    }
  mesh.level = median(mesh.mesh, mesh.nBlocksX * mesh.nBlocksY);
  if(isnan(mesh.level))
    {
      goto done;
    }
  interpolationAxis(n_x, mesh.nBlocksX, (int *)mesh.columnBlock, (double *)mesh.columnWeight);
  interpolationAxis(n_y, mesh.nBlocksY, (int *)mesh.rowBlock, (double *)mesh.rowWeight);
  /* an interpolation between blocks within half an ADU of their median rounds to 0 everywhere */
  for(k = 0; k < (long)mesh.nBlocksX * mesh.nBlocksY; k++)
    {
      if(fabs(mesh.mesh[k] - mesh.level) >= 0.5)
	{
	  overRanges(subtractMesh, &mesh, n_y, params->nThreads);
	  break;
	}
    }
  status = columnLevels(image, n_x, n_y, params->nThreads, columns, &levels->noise);
  if(status < 0)
    {
      *nSpots = (status == -1) ? SPOTS_BAD_IMAGE : SPOTS_NO_MEMORY;
      goto done;
    }
  levels->background = median(columns, n_x);
  if(isnan(levels->background))
    {
      goto done;
    }
  levels->threshold = floor(params->findSigma * levels->noise);
  if(!(levels->threshold >= MIN_THRESHOLD))
    {
      levels->threshold = MIN_THRESHOLD;
    }
  threshold = (int)levels->threshold;

  /* subtract the column levels, noting the pixels above the threshold */
  for(x = 0; x < n_x; x++)
    {
      rounded[x] = (int)floor(columns[x] + 0.5);
    }
  for(y = 0; y < n_y; y++)
    {
      int *row = image + (long)y * n_x;
      for(x = 0; x < n_x; x++)
	{
	  row[x] -= rounded[x];
	  if(row[x] > threshold)
	    {
	      if(nBright < MAX_BRIGHT_PIXELS)
		{
		  bright[nBright] = (long)y * n_x + x;
		}
	      nBright++;
	    }
	}
    }
  if(nBright > MAX_BRIGHT_PIXELS)
    {
      *nSpots = SPOTS_FLOODED;
      goto done;
    }

  /* group them into regions */
  for(k = 0; k < nBright && nRegions <= MAX_SPOTS; k++)
    {
      if(!mask[bright[k]])
	{
	  region *found = &regions[nRegions];
	  const int size = floodRegion(image, mask, n_x, n_y, bright[k], threshold, &stack, &capacity, found);
	  if(size < 0)
	    {
	      goto done;
	    }
	  if(size >= params->nmin && found->x > margin && found->x < n_x - margin
	     && found->y > margin && found->y < n_y - margin)
	    {
	      nRegions++;
	    }
	}
    }
  if(nRegions > MAX_SPOTS)
    {
      *nSpots = SPOTS_TOO_MANY;
      goto done;
    }

  context.image = image;
  context.n_x = n_x;
  context.n_y = n_y;
  context.params = params;
  context.regions = regions;
  overRanges(centroidRegions, &context, nRegions, params->nThreads);

  /* of regions that settled on the same spot, keep the one with the highest peak */
  dropped = calloc(nRegions > 0 ? nRegions : 1, sizeof(unsigned char));
  if(dropped == NULL)
    {
      goto done;
    }
  byY = malloc(sizeof(rankedY) * (nRegions > 0 ? nRegions : 1));
  rankOf = malloc(sizeof(int) * (nRegions > 0 ? nRegions : 1));
  near = malloc(sizeof(int) * (nRegions > 0 ? nRegions : 1));
  if(byY == NULL || rankOf == NULL || near == NULL)
    {
      goto done;
    }
  for(a = 0; a < nRegions; a++)
    {
      byY[a].y = regions[a].y;
      byY[a].index = a;
    }
  qsort(byY, nRegions, sizeof(rankedY), compareRankedY);
  for(a = 0; a < nRegions; a++)
    {
      rankOf[byY[a].index] = a;
    }
  for(a = 0; a < nRegions; a++)
    {
      /* the regions after a, in their order, within closest of it in y */
      int nNear = 0, rank;
      for(rank = rankOf[a] - 1; rank >= 0 && regions[a].y - byY[rank].y < closest; rank--)
	{
	  if(byY[rank].index > a)
	    {
	      near[nNear++] = byY[rank].index;
	    }
	}
      for(rank = rankOf[a] + 1; rank < nRegions && byY[rank].y - regions[a].y < closest; rank++)
	{
	  if(byY[rank].index > a)
	    {
	      near[nNear++] = byY[rank].index;
	    }
	}
      qsort(near, nNear, sizeof(int), compareInt);

      for(k = 0; k < nNear && !dropped[a]; k++)
	{
	  b = near[k];
	  if(!dropped[b] && hypot(regions[a].x - regions[b].x, regions[a].y - regions[b].y) < closest)
	    {
	      dropped[regions[b].peak > regions[a].peak ? a : b] = 1;
	    }
	}
      nKept += !dropped[a];
    }

  spots = malloc(sizeof(spot) * (nKept > 0 ? nKept : 1));
  if(spots == NULL)
    {
      goto done;
    }
  nKept = 0;
  for(a = 0; a < nRegions; a++)
    {
      if(!dropped[a])
	{
	  spot *out = &spots[nKept++];
	  int nAperture = 0;

	  out->x = regions[a].x;
	  out->y = regions[a].y;
	  out->x2 = regions[a].x2;
	  out->y2 = regions[a].y2;
	  out->xy = regions[a].xy;
	  out->peak = regions[a].peak;
	  out->background = rounded[columnOf(out->x, n_x)]
	    + meshDeviation(&mesh, columnOf(out->x, n_x), columnOf(out->y, n_y));
	  out->flux = NAN;
	  out->fluxerr = NAN;
	  if(params->aperture > 0)
	    {
	      out->flux = apertureSum(image, n_x, n_y, out->x, out->y, params->aperture, &nAperture);
	      if(nAperture > 0)
		{
		  /* noise of the pixels summed, and photon noise of the flux */
		  out->fluxerr = sqrt(nAperture * levels->noise * levels->noise
				      + (out->flux > 0 && params->gain > 0 ? out->flux / params->gain : 0.));
		}
	    }
	}
    }
  *nSpots = nKept;

 done:
  free(mesh.mesh);
  free((int *)mesh.columnBlock);
  free((double *)mesh.columnWeight);
  free((int *)mesh.rowBlock);
  free((double *)mesh.rowWeight);
  free(columns);
  free(rounded);
  free(bright);
  free(mask);
  free(stack);
  free(regions);
  free(dropped);
  free(byY);
  free(rankOf);
  free(near);
  return spots;
}
