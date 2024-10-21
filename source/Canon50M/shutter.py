import sys
import getopt
import time
import fcntl
import struct
import os

# Updated command numbers based on the C program output
GP0_OUTPUT_ENABLE = 0x40048100
GP0_SET_VALUE = 0x40048101
GP0_GET_VALUE = 0x80048102
S_CLOSE = 0
S_OPEN = 1

def print_usage_syntax(program_name):
    print(f"""Metrology camera shutter control.
Usage: {program_name} [options...]
    -o, --open     Open shutter.
    -c, --close    Close shutter.
    -t, --time     Set timing for shutter (ms).
    -d, --device   Specify device path (default: /dev/ttyUSB0).
    -h, --help     Display help message.
    -v, --verbose  Turn on verbose.""", file=sys.stderr)

def ioctl_command(fd, command, value):
    try:
        result = fcntl.ioctl(fd, command, struct.pack('I', value))
        return True
    except IOError as e:
        print(f"ioctl failed. Command: {hex(command)}, Value: {value}, Error: {e}")
        return False

def main(argv):
    verbose = False
    shutter = S_CLOSE
    timing = 0
    device_path = "/dev/ttyUSB0"

    try:
        opts, args = getopt.getopt(argv[1:], "oct:d:vh", ["open", "close", "time=", "device=", "verbose", "help"])
    except getopt.GetoptError:
        print_usage_syntax(argv[0])
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-o", "--open"):
            shutter = S_OPEN
        elif opt in ("-c", "--close"):
            shutter = S_CLOSE
        elif opt in ("-t", "--time"):
            timing = int(arg)
            shutter = -1
        elif opt in ("-d", "--device"):
            device_path = arg
        elif opt in ("-v", "--verbose"):
            verbose = True
        elif opt in ("-h", "--help"):
            print_usage_syntax(argv[0])
            sys.exit()

    if len(argv) < 2:
        print_usage_syntax(argv[0])
        sys.exit(2)

    if timing != 0 and timing < 10:
        print("Number too small. The unit of time is ms!!!")
        sys.exit(2)

    try:
        fd = os.open(device_path, os.O_RDWR)
        print(f"Successfully opened device: {device_path}")

        # Enable GPIO_1 to output mode
        if not ioctl_command(fd, GP0_OUTPUT_ENABLE, 1):
            raise IOError("Failed to set GPIO to output mode")

        if verbose:
            print(f"GP0_OUTPUT_ENABLE = {GP0_OUTPUT_ENABLE} {hex(GP0_OUTPUT_ENABLE)}")
            print(f"GP0_SET_VALUE = {GP0_SET_VALUE} {hex(GP0_SET_VALUE)}")
            print(f"GP0_GET_VALUE = {GP0_GET_VALUE} {hex(GP0_GET_VALUE)}")

        if shutter == S_OPEN:
            print("Opening shutter.")
            if not ioctl_command(fd, GP0_SET_VALUE, S_OPEN):
                raise IOError("Failed to open shutter")
        elif shutter == S_CLOSE:
            print("Closing shutter.")
            if not ioctl_command(fd, GP0_SET_VALUE, S_CLOSE):
                raise IOError("Failed to close shutter")
        elif timing != 0:
            print(f"Opening shutter for {timing} msec.")
            if not ioctl_command(fd, GP0_SET_VALUE, S_OPEN):
                raise IOError("Failed to open shutter")
            time.sleep(timing / 1000)
            if not ioctl_command(fd, GP0_SET_VALUE, S_CLOSE):
                raise IOError("Failed to close shutter")

    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)
    finally:
        os.close(fd)
        print(f"Closed device: {device_path}")

if __name__ == "__main__":
    main(sys.argv)