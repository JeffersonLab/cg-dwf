#
# The main header files. dwf-tables.h 
# I put dwf.c here because it contains macros and inline stuff and is really meand for
# inclusion. Perhaps dwf.h would have been a better name
CODE_HDRS = dwf-tables.h \
	    dwf.c

# --------------- SSE headers and sources ----------------
DWF_SSE_HDRS = dwf-ssef.h \
	       dwf-ssed.h

DWF_SSE_SRCS = dwf-ssef.c \
	       dwf-ssed.c
# --------------------------------------------------------

# -------------- BG/L Bluelight headers and sources ------
DWF_BLUELIGHT_HDRS = dwf-bluelightd.h \
	             dwf-bluelightf.h

DWF_BLUELIGHT_SRCS = dwf-bluelightf.c \
		     dwf-bluelightd.c
# --------------------------------------------------------

# -------------- Altivec headers and sources -------------
DWF_ALTIVEC_HDRS = dwf-altivecf.h
DWF_ALTIVEC_SRCS = dwf-altivecf.c
# --------------------------------------------------------

# Library name declaration
lib_LIBRARIES=libcg-dwf.a
libcg_dwf_a_SOURCES =

# Header file definition
include_HEADERS= $(CODE_HDRS)

# Architecture dependent construction of list of files to compile
#
# -----------------  SSE headers and sources --------------------
if ARCH_SSE
include_HEADERS  += $(DWF_SSE_HDRS)
libcg_dwf_a_SOURCES += $(DWF_SSE_SRCS)
endif
# ---------------------------------------------------------------

# -----------------  BG/L Bluelight headers and sources ---------
if ARCH_BLUELIGHT
include_HEADERS  += $(DWF_BLUELIGHT_HDRS)
libcg_dwf_a_SOURCES += $(DWF_BLUELIGHT_SRCS)
endif
# ---------------------------------------------------------------

# ----------------  AltiVec headers and sources -----------------
if ARCH_ALTIVEC
include_HEADERS  += $(DWF_ALTIVEC_HDRS)
libcg_dwf_a_SOURCES += $(DWF_ALTIVEC_SRCS)
endif
# ---------------------------------------------------------------

bin_SCRIPTS = cg-dwf-config
# Set up flags:

# General flags
INCFLAGS= -I@top_srcdir@

# Architecture dependent flags
DWF_CFLAGS=@DWF_CFLAGS@

# QMP flags
QMP_CFLAGS=@QMP_CFLAGS@

# Add it all onto the automake list of flags
AM_CFLAGS=$(INCFLAGS) $(DWF_CFLAGS) $(QMP_CFLAGS)
