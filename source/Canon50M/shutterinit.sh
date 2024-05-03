#! /bin/bash
rmmod pl2303
modprobe usbserial
insmod $ICS_MCSACTOR_DIR/source/ShutterDriver/pl2303.ko
usermod -a -G dialout pfs-data
shutter -o
shutter -c

