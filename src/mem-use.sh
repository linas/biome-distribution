#! /bin/sh

vmstat 60 >> mem.log &

while /bin/true ; do
	# echo the time is
	date >> mem.log
	ps aux |grep guile |grep -v grep >> mem.log
	sleep 300
done
