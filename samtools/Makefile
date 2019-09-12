CC=			gcc
CFLAGS=		-O2
DFLAGS=		-D_FILE_OFFSET_BITS=64 -D_USE_KNETFILE -D_CURSES_LIB=1
KNETFILE_O=	knetfile.o
LOBJS=		bgzf.o kstring.o bam_aux.o bam.o bam_import.o sam.o bam_index.o	\
			bam_pileup.o bam_lpileup.o bam_md.o glf.o razf.o faidx.o \
			$(KNETFILE_O) bam_sort.o sam_header.o bam_reheader.o kprobaln.o
AOBJS=		bam_tview.o bam_maqcns.o bam_plcmd.o sam_view.o	\
			bam_rmdup.o bam_rmdupse.o bam_mate.o bam_stat.o bam_color.o	\
			bamtk.o kaln.o bam2bcf.o bam2bcf_indel.o errmod.o sample.o
PROG=		samtools
INCLUDES=	-I.
SUBDIRS=	.
LIBPATH=
LIBCURSES=	-lcurses # -lXCurses

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all-recur lib-recur clean-recur cleanlocal-recur install-recur:
		@target=`echo $@ | sed s/-recur//`; \
		wdir=`pwd`; \
		list='$(SUBDIRS)'; for subdir in $$list; do \
			cd $$subdir; \
			$(MAKE) CC="$(CC)" DFLAGS="$(DFLAGS)" CFLAGS="$(CFLAGS)" \
				INCLUDES="$(INCLUDES)" LIBPATH="$(LIBPATH)" $$target || exit 1; \
			cd $$wdir; \
		done;

all:libbam.a

.PHONY:all lib clean cleanlocal
.PHONY:all-recur lib-recur clean-recur cleanlocal-recur install-recur

lib:libbam.a

libbam.a:$(LOBJS)
		$(AR) -cru $@ $(LOBJS)

cleanlocal:
		rm -fr gmon.out *.o a.out *.exe *.dSYM razip bgzip $(PROG) *~ *.a *.so.* *.so *.dylib

clean:cleanlocal-recur
