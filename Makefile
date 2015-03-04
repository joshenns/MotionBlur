all: motionblur

motionblur: motionblur.cc fft.cc
	g++ -I /usr/include/netpbm -o motionblur motionblur.cc fft.cc -Wall -lnetpbm

clean:
	rm -f *~

clean-all: clean
	rm -f motionblur
