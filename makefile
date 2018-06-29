CXX=g++ -O2 -Wall
### -O2 -O5 -Os
#g++ `pkg-config --cflags itpp` -o hello.out hello.cpp `pkg-config --libs itpp`

START=`pkg-config --cflags itpp`
END=`pkg-config --libs itpp`
files=mm_read.c mm_read.h mmio.c mmio.h mm_write.c mm_write.h lib.cpp lib.h my_lib.h makefile
command=$(CXX) $(START) -o $@ $< $(word 2,$^) $(word 4, $^) $(word 6, $^) $(word 8, $^) $(END)

all: parserTest.out ldpcTest.out test.out alist_test.out
###include all headfiles into my_lib.h
bp_decoding3.out:bp_decoding3.c mm_read.c mm_read.h mmio.c mmio.h mm_write.c mm_write.h lib.cpp lib.h my_lib.h makefile
	$(command)
data_collect_perm_check.out:data_collect_perm_check.c $(files)
	$(command)
convert_data.out:convert_data.c mm_read.c mm_read.h mmio.c mmio.h mm_write.c mm_write.h my_lib.h
	$(CXX) $(START) -o $@ $< $(word 2,$^) $(word 4, $^) $(word 6, $^) $(END)
error_analysis.out:error_analysis.c mm_read.c mm_read.h mmio.c mmio.h mm_write.c mm_write.h my_lib.h
	$(CXX) $(START) -o $@ $< $(word 2,$^) $(word 4, $^) $(word 6, $^) $(END)
rand_decode_perm_check.out:rand_decode_perm_check.c mm_read.c mm_read.h mmio.c mmio.h mm_write.c mm_write.h lib.cpp lib.h my_lib.h
	$(CXX) $(START) -o $@ $< $(word 2,$^) $(word 4, $^) $(word 6, $^) $(word 8, $^) $(END)
rand_decode_perm.out:rand_decode_perm.c mm_read.c mm_read.h mmio.c mmio.h mm_write.c mm_write.h lib.cpp lib.h my_lib.h
	$(CXX) $(START) -o $@ $< $(word 2,$^) $(word 4, $^) $(word 6, $^) $(word 8, $^) $(END)
rand_decode3.out:rand_decode3.c mm_read.c mm_read.h mmio.c mmio.h mm_write.c mm_write.h lib.cpp lib.h my_lib.h
	$(CXX) $(START) -o $@ $< $(word 2,$^) $(word 4, $^) $(word 6, $^) $(word 8, $^) $(END)
rand_decode.out:rand_decode.c mm_read.c mm_read.h mmio.c mmio.h mm_write.c mm_write.h lib.cpp lib.h my_lib.h
	$(CXX) $(START) -o $@ $< $(word 2,$^) $(word 4, $^) $(word 6, $^) $(word 8, $^) $(END)
gauge_to_stabilizer.out:gauge_to_stabilizer.c mm_read.c mm_read.h mmio.c mmio.h mm_write.c mm_write.h my_lib.h
	$(CXX) $(START) -o $@ $< $(word 2,$^) $(word 4, $^) $(word 6, $^) $(END)
decoding.out:decoding.c mm_read.c mm_read.h mmio.c mmio.h mm_write.c mm_write.h my_lib.h
	$(CXX) $(START) -o $@ $< $(word 2,$^) $(word 4, $^) $(word 6, $^) $(END)
test.out:test.c mm_read.c mm_read.h mmio.c mmio.h mm_write.c mm_write.h lib.cpp lib.h my_lib.h
	$(CXX) $(START) -o $@ $< $(word 2,$^) $(word 4, $^) $(word 6, $^) $(word 8, $^) $(END)
code_generator.out:code_generator.c mmio.c mmio.h mm_write.c mm_write.h
	$(CXX) $(START) -o $@ $< $(word 2, $^) $(word 4, $^) $(END)
clean:
	rm  *~
	rm \#*

