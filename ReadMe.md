Using the method of four Russians for matrix multiplication in GF2, GF5, and GF7
with HPC considerations

I've included every file I bothered to keep this summer, so many are not usefull and or
not a full working implementation.

First here is a list of files containing unique working solutions, in order that I made them.
* full inline documentation

GF2
   fourRussian.c*
   splitN.c
   scriptTest3.c
   kernelt4, 8.c

GF5
   t2fourrusgf5.c (bit sliced)
   unrollt2, 4, 8.c (bit packed)
   bitslicedforust2.c* (bit sliced)
   paddedslicedintrinsicsT2.c (bit sliced)
   unrollslicet2, 4, 8.c (bit sliced)
   paralleAttempt2.c (bit sliced)

GF7 (all bit sliced)
   slicedt2gf7.c
   unrollintinsicT1, 2, 4.c
   endWithMod7T1, 2, 4.c

Full list of files 

GF2
   assembly.c
   assembly2.c
   b.c
   classical.c
   copy.c
   copy2.c
   copy3.c
   copy4.c
   copy5.c
   copy6.c
   copy7.c
   copy8.c
   fourRussian.c
   hopefullyTheLast.c
   kernal.c
   kernalt4.c
   kernalt8.c
   kernalt8n512.c
   kernalt8reverse.c
   main.c
   notes.txt
   optimization.c
   reassembly.c
   reorganization.c
   sampleRun.txt
   scriptTest.c
   scriptTest.sh
   scriptTest2.c
   scriptTest3.c
   speedTest.c
   splitN.c
   test.c
   timings.xlsx

GF5
   addition.c
   bitpackedAdditionGF5.c
   bitpackedAdditionGF7.c
   bitslicedforrust2.c
   checkGF7.c
   columnmajorApackedunrolled.c
   fourrusgf5.c
   kernalt2.c
   lunchScript.sh
   negation.c
   optimizeKernel.c
   packedbitsNegation.c
   paddedSlicedIntrinsicsT2.c
   parallelt2.c
   results.xlsx
   scrach.txt
   slicedintrinsicsT2.c
   stanfordMod5.c
   t1fourrusgf5.c
   t2fourrusgf5.c
   test.c
   timeApproximations.c
   unrollksliceT2.c
   unrollsliceT1.c
   unrollsliceT2.c
   unrollsliceT4.c
   unrollsliceTest.c
   unrollt1.c
   unrollt2.c
   unrollt4.c

GF7
   endWithMod7T1.c
   endWithMod7T2.c
   endWithMod7T4.c
   slicedt2gf7.c
   unrollintinsicT1.c
   unrollintinsicT2.c
   unrollintinsicT4.c

Here is a short description of each file, also included at the top of each file



GF2


assembly.c uses assembly macros for the
* loading and bitwise operations used in
* the matrix multiplication

assembly2.c is an older version of assebly.c
* except the createTable function doesn't compile.
* But there might be something of use here so I
* haven't deleted it

b.c was used for random testing to put later
* into other files, each test typically stands
* on its own and is commented out as a block.
* Some of these tests include: updating the creatTable
* function,
* properly padding matrixes, A, B, C

classical.c was my starting point and just initializes
* three matrixies and multiplies them without modulus
* then prints them out

copy.c has the first implementation of the method of four
* Russians. It seems to be limited to multiples of 32 for sizes
* of n, but accepts any value for t which was later taken out
* for simplicity.

copy2.c adds a testing method at the end of the file,
* but restricts t to a divisor of 32

copy3.c changes A to be stored in column major order
* A stores 32 elements from one row in an int then works
* down the rows 32 columns at a time

copy4.c adds subTables to the four Russian lookup table.
* Previously createTables only created tables for t rows at
* a time. Now it creates tables for 32 rows at a time since
* that many elemnts of A are read at one. subTables are 
* created t rows at a time and there are 32/t subTables.
*
* The outer most loop is now la (number of integer columns)
* Perviously it was tloops = lb/t (number of columns over t)

copy5.c unrolls the main loop doing the multiplication
* It doesn't work at sizes other than those set and so
* is a test and not a final product that should be used
* in any setting

copy6.c unrolls the matrix multiply loop further and
* uses instrinsic intructions. createTable use makes
* use of instrinsic intructions.
*
* This implementation doesn't correctly compute matrixes
* with dimensions, l, or n, not a multiple of 256

copy7.c adds extra loops around copy6.c to handle
* a matrix of any size by using 64 bit instructions
* on any overflow past multiples of 256 elementsv

copy8.c rewrites several instinsic instructions as
* assembly macros that compile better.
*
* The major change of this file is that is only computes
* the result for a matrixes where m, l, n = 4, 64, 256 
* and so is a just a µkernel
*
* For an unkown reason it runs ~50% of peak (16 cycles)
* the first time it's run and ~10% at other times

fourRussina.c is my first try at implementing the
* method of four Russians on GF2. It is the slowest
* and most general result I have since it accepts any
* value for t, but doesn't consider any HPC optimizations
*
* Although this program is one of my slowest most of the 
* concepts outside of optimization are best shown here,
* so I've included full documentation throughout this file.

hopefullyTheLast.c was my first attempt at a kernel 
* version of my program which I was hoping was the last,
* but wasn't.
*
* This implementation hits about 5% of peak (16 cycles).
* Later attempts use different unrolling and ways of 
* parsing indexs from A

kernal.c resembles my final immplementation, but
* runs with a bus error, so its probably not useful

kernalt4.c is one of my final kernals with t = 4
* it has similar performance to kernalt8.c and kernalt8reverse.c

kernalt8.c is one of my final kernals with t = 8
* it has similar performance to kernalt4.c and kernalt8reverse.c

kernalt8n512 is another final implementation of my kernal,
* with t = 8 and n = 512. Doubling n however didn't lead to 
* any higher performance instead lowering it.

kernalt4.c is one of my final kernals with t = 8
* it has similar performance to kernalt4.c and kernalt8.c
*
* This kernel differnts from the other by reading through
* integers of A and subTables in T in reverse.  This makes
* things run slightly faster, but more so if the compiler
* followed the code closer.

main.c immplements an early version of createTable, but doesn't use to to do a matrix multiplication

optimization.c attempts to use cache effectivly and unroll
* intrinsic instructions enough to saturate the pipeline.
* This makes it one of my best implementations for all matrix
* sizes in GF2

reassembly.c rewrites intrinsic instructions as macros which are compiled better
* it runs with a seg fault so don't use this file

reorganization.c reorganizies the loops of the matrix multiply
* for better cache use.  The order is over integer columns of A,
* rows of A and then integer columns of B.

scriptTest 1, 2, 3.c are versions that compute after taking in
* each dimension as a command line argument.
* They work with scriptTest.sh to test what sizes they correclty
* compute
*
* run by: ./a.out 512 512 512

speedTest.c is an early immplementation I run speed tests on

splitN.c splits up computation over 256 elements of N
* and moves the n loop to the outer most.
* This limits how large the lookup table can get and shouldn't
* harm performance for small sizes, but does slightly.

GF2/test.c contants some early tests, but eventually I made my fastest full implementation here
* that isn't a kernel




GF5


addition.c was my first attempt at addition and modulus
* in GF5 using bit sliced element storage

bitpackedAdditionGF7.c was my first attempt for addition and modulus
* for bitpacked GF5 elements. It doesn't work but was fixed in full
* matrix multiply programs.

bitpackedAdditionGF7.c was my first attempt for addition and modulus
* for bitpacked GF7 elements. It doesn't work and was never fixed
* since a bitpacked four Russian matrix multiply was never made.

bitslicedforrust2.c is a bitsliced Four Russian matrix multiply with t = 2
* Since it isn't unrolled and doesn't use intrinsics I've included full 
* documentation in this file

checkGF7.c is a addition and modulus of bit sliced
* elements of GF7

columnmajorApackedunrolled.c is equivalent to unrollslice.c files
* except A is store in column major order. It hasn't been amended from 
* a row major A except to see if reading from a column major A helps.
* It doesn't help so I never made a full columnmajor A implementation.

fourrusgf5.c was my first attempt at a four Russina matrix multiply
* in GF5. It uses bitpacked matrixes and so doesn't run as fast as sliced
* for sizes over 64.

kernalt2.c is a kernel with m, l, n sizes 2, 16, 64
* It reaches ~20% of peak (22 cycles), but this is for
* a bitpacked operation and slicing would be ~2.5x faster
* A bitsliced kernel was never made

negation.c has a few different operations for GF5 commented out
* examples: doubling, negation, substraction, addition

optimizeKernel.c is an extension of kernalt2.c,
* but doens't run much quicker
*
* Again, this is for a bitpacked operation and slicing
* would be ~2.5x faster. A bitsliced kernel was never made

packedbitsNegation.c gives bitwise negations and doubling for packed
* elements of GF5

paddedSlicedIntrinsicsT2.c added padding to the first intrinsic bitslcied
* implementation so it can successfully compute all sizes of matrixes

parallelt2.c First attemp at adding parallism with OpenMP
* Doesn't do anything

parallelAttempt.c add a parallel for loop added to the loop going through rows of A
* Either the OpenMP overhead or cache conflit results in 16 threads (on a 32 threads CPU)
* having the same performance as one. 32 and 4 threads have sever times worse performance

slicedintrinsicsT2.c changed createTable and the main loop to
* instrinsic instructions

stanfordMod5.c taken from stanford's bittwiddling hacks website
* not used in further work since it doesn't scale to registers
* with multiple elements as far as I can tell.

t1fourrusgf5.c is the same implementation as fourrusgf5.c
* but with t = 1

t2fourrusgf5.c is the same implementation as fourrusgf5.c
* but with t = 2

GF5/test.c is just a sample row computation I did to varify my testing routine

timeAproximations.c gives an order of magnitude approximation of
* how long each field needs for a four Russian matrix multiplication

unrollksliceT2.c tries to speedup by unrolling the loop
* that goes through the different subTables (whose int is k)
* It doesn't work as well as the other loop I tried unrolling

unrollsliceT 1, 2, 4.c are my fastest GF5 implementation with
* t = 1, 2, 4

unrollsliceTest.c tested whether unrolling the loop twice acts
* better than a 4x unroll.
* Result: they take the same amount of time
* So the pipeline is probably saturated with two and four times unroll

unrollt 1, 2, 4.c are unrolled bitpacked four Russian matrix multiplies.
* Since the compiler really likes t*fourrusgf5.c and uses
* instrinsics and unrolls automatically, unrollt*.c runs slower.





GF7


endWithMod7T1, 2, 4.c are exactly unrollintinsicT*.c except
* the elements in C are allowed to be 7 throughout the computation
* and then C is read through again and 7's are changed to 0's
* (This speeds up the program)

slicedt2gf7.c is a modification of t2fourrusgf5.c to work
* with GF7 addition

unrollintinsicT 1, 2, 4.c are modifications of unrollt*.c to work
* with GF7 additions


























