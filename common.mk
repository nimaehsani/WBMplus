UNIX=$(shell uname)
ifndef GHAASDIR
export GHAASDIR=/usr/local/share/ghaas
endif

export UNIXCC=gcc
export UNIXCCOPS=-g -Wall -pthread -fsigned-char -D_GNU_SOURCE
export UNIXLIBS=-pthread -lm
export UNIXMAKE=make