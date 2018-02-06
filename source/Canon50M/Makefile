# Comment/uncomment the following line to disable/enable debugging
DEBUG = y

INSTALLPATH = /usr/local/bin

VER_SERIES=$(shell uname -r | cut -f1-2 -d.)
VER_MINOR=$(shell uname -r | cut -f3 -d.)

CC=gcc

DEBFLAGS = -fPIC -O2 -Wno-format-zero-length -Wno-deprecated-declarations

ALL_CFLAGS += $(DEBFLAGS) $(DO_MLOCK) $(BUILDOPTIONS)
LDFLAGS += $(ALL_CFLAGS) -L. -L/opt/EDTpdv/ -L/usr/local/lib  -I/cfht/include -I. -I/usr/local/include  -I/opt/EDTpdv/
LDLIBS += -lpdv -lpthread -lm -ldl -lcfitsio


TARGET =  canonexp shutter canoninit


%.o : %.c
	@echo $(CC) $(ALL_CFLAGS) -c $< $(LOGOUT)
	@$(CC) $(ALL_CFLAGS) -c $< $(LOGOUT)
	@$(SHOWOUT)


all:  $(TARGET)

	
install: $(TARGET)
	@for file in $(TARGET); do \
      echo "--> Installing scripts $$file";\
      sudo cp  -fp $$file $(INSTALLPATH)/$$file-$(shell date +%y%m%d) ;\
      sudo cp  -fp $$file $(INSTALLPATH) ;\
	done

clean:
	rm -rf *.o *.fits
	@for file in $(TARGET); do \
	  rm -rf $$file ;\
	done

