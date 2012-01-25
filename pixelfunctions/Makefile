#!/usr/bin/make -f

.PHONY: all clean check dist

OBJS = pixfunplugin.o pixelfunctions.o
CFLAGS := -fPIC -Wall -Wno-long-long -pedantic \
          $(shell gdal-config --cflags) $(CFLAGS)

#CFLAGS := -O0 -g3 -ggdb3 $(CFLAGS)
CFLAGS := -O3  $(CFLAGS)

ARCHIVE = pixfun-plugin-$(shell date +%Y%m%d)

TARGET = gdal_PIXFUN.so
#TARGET = gdal_PIXFUN.dylib

all: $(TARGET)

clean:
	$(RM) $(TARGET) *.o *~

dist:
	$(RM) $(ARCHIVE).tar.gz
	mkdir -p $(ARCHIVE)/tests/data
	cp $(OBJS:.o=.c) Makefile README.txt $(ARCHIVE)
	cp tests/*.py $(ARCHIVE)/tests
	cp tests/data/*.vrt tests/data/*.tif $(ARCHIVE)/tests/data
	tar cvfz $(ARCHIVE).tar.gz $(ARCHIVE)
	$(RM) -r $(ARCHIVE)

check: $(TARGET)
	cd tests && python test_pixfun.py

$(TARGET): $(OBJS)
	$(CC) -shared -o $@ $(OBJS) $(shell gdal-config --libs)
