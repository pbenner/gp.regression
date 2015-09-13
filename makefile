# Makefile for the R library
R          ?= R
RFLAGS     ?= --no-multiarch
VERSION    ?= 1.0-1
PACKAGE    ?= gp.regression
TARGET      = $(PACKAGE)_$(VERSION).tar.gz

all: $(TARGET)

$(TARGET): NAMESPACE
	$(RM) .Rhistory
	$(RM) */.Rhistory
	(cd .. && $(R) CMD build $(RFLAGS) $(PACKAGE))

.SECONDEXPANSION:
NAMESPACE: $$(wildcard R/*.R) $$(wildcard src/*.c)
	echo "library(roxygen2); roxygenize('.')" | $(R) --vanilla --quiet
	touch $@

check: $(TARGET)
	$(R) CMD check $(RFLAGS) ../$(TARGET)

install: $(TARGET)
	$(R) CMD INSTALL $(RFLAGS) ../$(TARGET)

clean:
	$(RM) -r $(PACKAGE).Rcheck
	$(RM) -r $(PACKAGE).pdf
	$(RM) -r man
	$(RM) ../$(TARGET)
	$(RM) NAMESPACE
