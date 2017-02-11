#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
# Build supernova
#

LIBPY=lib/python
ASSEMBLYBINS=CP DF ParseBarcodedFastqs FastFastbCount MakeFasta
ASSEMBLYPATH=$(shell pwd)/lib/assembly
VERSION=$(shell git describe --tags --always --dirty)
JEMALLOC_ROOT := $(shell pwd)/lib/jemalloc
JEMALLOC_SRC := $(JEMALLOC_ROOT)/build
JEMALLOC_CONF := $(JEMALLOC_SRC)/configure
export JEMALLOC_VERSION=3.6.0
export JEMALLOC_PATH := $(JEMALLOC_ROOT)/$(JEMALLOC_VERSION)
export DEPEND_FLAGS := -static
export GOPATH=$(shell pwd)/lib/go
export CS =

.PHONY: $(ASSEMBLYBINS) test clean jemalloc

#
# Targets for development builds.
#
all: $(JEMALLOC_PATH) $(ASSEMBLYBINS) 

$(JEMALLOC_PATH): $(JEMALLOC_CONF)
	make -C $(JEMALLOC_SRC)
	make -C $(JEMALLOC_SRC) install_lib install_bin install_include DESTDIR=..

$(JEMALLOC_CONF):
	cd $(JEMALLOC_SRC); \
	./autogen.sh ;\
	CPPFLAGS="-Wno-unused -Wno-maybe-uninitialized" ./configure --prefix=/$(JEMALLOC_VERSION) --with-jemalloc-prefix=je_

jemalloc: $(JEMALLOC_PATH)

lib/bin:
	mkdir lib/bin

$(ASSEMBLYBINS): lib/bin $(ASSEMBLYPATH)
	$(MAKE) -j16 -C $(ASSEMBLYPATH) $@ STATIC=yes
	cp $(ASSEMBLYPATH)/bin/$@ lib/bin


test:
	MROPATH=$(PWD)/mro \
	    PATH=$(PATH):$(PWD)/bin \
	    PYTHONPATH=$(PYTHONPATH):$(PWD)/lib/python \
	    nosetests --with-xunit --with-coverage --cover-erase --cover-html --cover-xml --cover-package=stages,tenkit,kitten mro/stages/* lib/python

clean: 
	rm -rf $(GOPATH)/bin
	rm -rf $(GOPATH)/pkg
	rm -rf $(ASSEMBLYPATH)/bin
	rm -rf lib/bin
	rm -rf $(JEMALLOC_PATH)
	rm -rf $(JEMALLOC_CONF)
	make -C lib/assembly cleanAll

