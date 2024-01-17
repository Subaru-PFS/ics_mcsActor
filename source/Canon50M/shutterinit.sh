#! /bin/bash
rmmod pl2303
modprobe usbserial
insmod /home/pfs/mhs/devel/ics_mcsActor/source/ShutterDriver/pl2303.ko
usermod -a -G dialout pfs-data
shutter -o
shutter -c

