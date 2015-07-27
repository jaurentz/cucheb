!-----------------------------------------------------------------------
      subroutine coocsr(nrow,nnz,a,ir,jc,ao,jao,iao)
!-----------------------------------------------------------------------
      real*8 a(*),ao(*),x
      integer ir(*),jc(*),jao(*),iao(*)
!-----------------------------------------------------------------------
!  Coordinate     to   Compressed Sparse Row
!-----------------------------------------------------------------------
! converts a matrix that is stored in coordinate format
!  a, ir, jc into a row general sparse ao, jao, iao format.
!
! on entry:
!---------
! nrow      = dimension of the matrix
! nnz      = number of nonzero elements in matrix
! a,
! ir,
! jc    = matrix in coordinate format. a(k), ir(k), jc(k) store the nnz
!         nonzero elements of the matrix with a(k) = actual real value of
!         the elements, ir(k) = its row number and jc(k) = its column
!        number. The order of the elements is arbitrary.
!
! on return:
!-----------
! ir       is destroyed
!
! ao, jao, iao = matrix in general sparse matrix format with ao
!       continung the real values, jao containing the column indices,
!      and iao being the pointer to the beginning of the row,
!      in arrays ao, jao.
!
! Notes:
!------ This routine is NOT in place.  See coicsr
!
!------------------------------------------------------------------------
      do 1 k=1,nrow+1
         iao(k) = 0
 1    continue
! determine row-lengths.
      do 2 k=1, nnz
         iao(ir(k)) = iao(ir(k))+1
 2    continue
! starting position of each row..
      k = 1
      do 3 j=1,nrow+1
         k0 = iao(j)
         iao(j) = k
         k = k+k0
 3    continue
! go through the structure  once more. Fill in output matrix.
      do 4 k=1, nnz
         i = ir(k)
         j = jc(k)
         x = a(k)
         iad = iao(i)
         ao(iad) =  x
         jao(iad) = j
         iao(i) = iad+1
 4    continue
! shift back iao
      do 5 j=nrow,1,-1
         iao(j+1) = iao(j)
 5    continue
      iao(1) = 1
      return
!------------- end of coocsr -------------------------------------------
!-----------------------------------------------------------------------
      end

      subroutine dcsort(ival, n, icnt, index, ilo, ihi)
!-----------------------------------------------------------------------
!     Specifications for arguments:
!     ----------------------------
      integer n, ilo, ihi, ival(n), icnt(ilo:ihi), index(n)
!-----------------------------------------------------------------------
!    This routine computes a permutation which, when applied to the
!    input vector ival, sorts the integers in ival in descending
!    order.  The permutation is represented by the vector index.  The
!    permuted ival can be interpreted as follows:
!      ival(index(i-1)) .ge. ival(index(i)) .ge. ival(index(i+1))
!
!    A specialized sort, the distribution counting sort, is used 
!    which takes advantage of the knowledge that
!        1)  The values are in the (small) range [ ilo, ihi ]
!        2)  Values are likely to be repeated often
!
!    contributed to SPARSKIT by Mike Heroux. (Cray Research) 
!    --------------------------------------- 
!----------------------------------------------------------------------- 
! Usage:
!------ 
!     call dcsort( ival, n, icnt, index, ilo, ihi )
!
! Arguments:
!----------- 
!    ival  integer array (input)
!          On entry, ia is an n dimensional array that contains
!          the values to be sorted.  ival is unchanged on exit.
!
!    n     integer (input)
!          On entry, n is the number of elements in ival and index.
!
!    icnt  integer (work)
!          On entry, is an integer work vector of length 
!          (ihi - ilo + 1).
!
!    index integer array (output)
!          On exit, index is an n-length integer vector containing
!          the permutation which sorts the vector ival.
!
!    ilo   integer (input)
!          On entry, ilo is .le. to the minimum value in ival.
!
!    ihi   integer (input)
!          On entry, ihi is .ge. to the maximum value in ival.
!
! Remarks:
!--------- 
!         The permutation is NOT applied to the vector ival.
!
!----------------------------------------------------------------
!
! Local variables:
!    Other integer values are temporary indices.
!
! Author: 
!-------- 
!    Michael Heroux
!    Sandra Carney
!       Mathematical Software Research Group
!       Cray Research, Inc.
!
! References:
!    Knuth, Donald E., "The Art of Computer Programming, Volume 3:
!    Sorting and Searching," Addison-Wesley, Reading, Massachusetts,
!    1973, pp. 78-79.
!
! Revision history:
!    05/09/90: Original implementation.  A variation of the 
!              Distribution Counting Sort recommended by
!              Sandra Carney. (Mike Heroux)
!
!-----------------------------------------------------------------
!     ----------------------------------
!     Specifications for local variables
!     ----------------------------------
      integer i, j, ivalj
!
!     --------------------------
!     First executable statement
!     --------------------------
      do 10 i = ilo, ihi
        icnt(i) = 0
 10   continue
!
      do 20 i = 1, n
        icnt(ival(i)) = icnt(ival(i)) + 1
 20   continue
!
      do 30 i = ihi-1,ilo,-1
        icnt(i) = icnt(i) + icnt(i+1)
 30   continue
!
      do 40 j = n, 1, -1
        ivalj = ival(j)
        index(icnt(ivalj)) = j
        icnt(ivalj) = icnt(ivalj) - 1
 40   continue
      return
      end
!-------end-of-dcsort---------------------------------------------------

!-----------------------------------------------------------------------
      subroutine csrjad(nrow, a, ja, ia, idiag, iperm, ao, jao, iao) 
      integer ja(*), jao(*), ia(nrow+1), iperm(nrow), iao(nrow) 
      real*8 a(*), ao(*)
!-----------------------------------------------------------------------
!    Compressed Sparse Row  to   Jagged Diagonal storage. 
!----------------------------------------------------------------------- 
! this subroutine converts  matrix stored in the compressed sparse
! row format to the jagged diagonal format. The data structure
! for the JAD (Jagged Diagonal storage) is as follows. The rows of 
! the matrix are (implicitly) permuted so that their lengths are in
! decreasing order. The real entries ao(*) and their column indices 
! jao(*) are stored in succession. The number of such diagonals is idiag.
! the lengths of each of these diagonals is stored in iao(*).
! For more details see [E. Anderson and Y. Saad,
! ``Solving sparse triangular systems on parallel computers'' in
! Inter. J. of High Speed Computing, Vol 1, pp. 73-96 (1989).]
! or  [Y. Saad, ``Krylov Subspace Methods on Supercomputers''
! SIAM J. on  Stat. Scient. Comput., volume 10, pp. 1200-1232 (1989).]
!----------------------------------------------------------------------- 
! on entry:
!---------- 
! nrow 	  = row dimension of the matrix A.
!
! a, 
! ia, 
! ja      = input matrix in compressed sparse row format. 
!
! on return: 
!----------
! 
! idiag = integer. The number of jagged diagonals in the matrix.
!
! iperm = integer array of length nrow containing the permutation
!         of the rows that leads to a decreasing order of the
!         number of nonzero elements.
!
! ao    = real array containing the values of the matrix A in 
!         jagged diagonal storage. The j-diagonals are stored
!         in ao in sequence. 
!
! jao   = integer array containing the column indices of the 
!         entries in ao.
!
! iao   = integer array containing pointers to the beginning 
!         of each j-diagonal in ao, jao. iao is also used as 
!         a work array and it should be of length n at least.
!
!----------------------------------------------------------------------- 
!     ---- define initial iperm and get lengths of each row
!     ---- jao is used a work vector to store tehse lengths
!     
      idiag = 0
      ilo = nrow 
      do 10 j=1, nrow
         iperm(j) = j 
         len = ia(j+1) - ia(j)
         ilo = min(ilo,len) 
         idiag = max(idiag,len) 
         jao(j) = len
 10   continue 
!     
!     call sorter to get permutation. use iao as work array.
!    
      call dcsort (jao, nrow, iao, iperm, ilo, idiag) 
!     
!     define output data structure. first lengths of j-diagonals
!     
      do 20 j=1, nrow
         iao(j) = 0
 20   continue
      do 40 k=1, nrow
         len = jao(iperm(k)) 
         do 30 i=1,len
            iao(i) = iao(i)+1
 30      continue
 40   continue
!     
!     get the output matrix itself
!     
      k1 = 1
      k0 = k1
      do 60 jj=1, idiag
         len = iao(jj)
         do 50 k=1,len
            i = ia(iperm(k))+jj-1
            ao(k1) = a(i)
            jao(k1) = ja(i) 
            k1 = k1+1
 50      continue
         iao(jj) = k0
         k0 = k1
 60   continue
      iao(idiag+1) = k1
      return
!----------end-of-csrjad------------------------------------------------
!-----------------------------------------------------------------------
      end

      subroutine aplb(nrow,ncol,job,a,ja,ia,b,jb,ib,c,jc,ic,nzmax,iw,
     *ierr)
      real*8 a(*), b(*), c(*) 
      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(nrow+1),ic(nrow+1),
     *iw(ncol)
!-----------------------------------------------------------------------
! performs the matrix sum  C = A+B. 
!-----------------------------------------------------------------------
! on entry:
! ---------
! nrow	= integer. The row dimension of A and B
! ncol  = integer. The column dimension of A and B.
! job   = integer. Job indicator. When job = 0, only the structure
!                  (i.e. the arrays jc, ic) is computed and the
!                  real values are ignored.
!
! a,
! ja,
! ia   = Matrix A in compressed sparse row format.
! 
! b, 
! jb, 
! ib	=  Matrix B in compressed sparse row format.
!
! nzmax	= integer. The  length of the arrays c and jc.
!         amub will stop if the result matrix C  has a number 
!         of elements that exceeds exceeds nzmax. See ierr.
! 
! on return:
!----------
! c, 
! jc, 
! ic	= resulting matrix C in compressed sparse row sparse format.
!	    
! ierr	= integer. serving as error message. 
!         ierr = 0 means normal return,
!         ierr .gt. 0 means that amub stopped while computing the
!         i-th row  of C with i=ierr, because the number 
!         of elements in C exceeds nzmax.
!
! work arrays:
!------------
! iw	= integer work array of length equal to the number of
!         columns in A.
!
!-----------------------------------------------------------------------
      logical values
      values = (job .ne. 0) 
      ierr = 0
      len = 0
      ic(1) = 1 
      do 1 j=1, ncol
         iw(j) = 0
 1    continue
!     
      do 500 ii=1, nrow
!     row i 
         do 200 ka=ia(ii), ia(ii+1)-1 
            len = len+1
            jcol    = ja(ka)
            if (len .gt. nzmax) goto 999
            jc(len) = jcol 
            if (values) c(len)  = a(ka) 
            iw(jcol)= len
 200     continue
!     
         do 300 kb=ib(ii),ib(ii+1)-1
            jcol = jb(kb)
            jpos = iw(jcol)
            if (jpos .eq. 0) then
               len = len+1
               if (len .gt. nzmax) goto 999
               jc(len) = jcol
               if (values) c(len)  = b(kb)
               iw(jcol)= len
            else
               if (values) c(jpos) = c(jpos) + b(kb)
            endif
 300     continue
         do 301 k=ic(ii), len
           iw(jc(k)) = 0
 301     continue
         ic(ii+1) = len+1
 500  continue
      return
 999  ierr = ii
      return
!------------end of aplb ----------------------------------------------- 
!-----------------------------------------------------------------------
      end

!----------------------------------------------------------------------- 
      subroutine vperm(n, x, perm) 
      integer n, perm(n) 
      real*8 x(n)
!-----------------------------------------------------------------------
! this subroutine performs an in-place permutation of a real vector x 
! according to the permutation array perm(*), i.e., on return, 
! the vector x satisfies,
!
!	x(perm(j)) :== x(j), j=1,2,.., n
!
!-----------------------------------------------------------------------
! on entry:
!---------
! n 	= length of vector x.
! perm 	= integer array of length n containing the permutation  array.
! x	= input vector
!
! on return:
!---------- 
! x	= vector x permuted according to x(perm(*)) :=  x(*)
!
!----------------------------------------------------------------------c
!           Y. Saad, Sep. 21 1989                                      c
!----------------------------------------------------------------------c
! local variables 
      real*8 tmp, tmp1
!
      init      = 1
      tmp       = x(init)
      ii        = perm(init)
      perm(init)= -perm(init)
      k         = 0
!     
! loop
! 
 6    k = k+1
!
! save the chased element --
! 
      tmp1      = x(ii) 
      x(ii)     = tmp
      next      = perm(ii) 
      if (next .lt. 0 ) goto 65
!     
! test for end 
!
      if (k .gt. n) goto 101
      tmp       = tmp1
      perm(ii)  = - perm(ii)
      ii        = next 
!
! end loop 
!
      goto 6
!
! reinitilaize cycle --
!
 65   init      = init+1
      if (init .gt. n) goto 101
      if (perm(init) .lt. 0) goto 65
      tmp       = x(init)
      ii        = perm(init)
      perm(init)=-perm(init)
      goto 6
!     
 101  continue
      do 200 j=1, n
         perm(j) = -perm(j)
 200  continue 
!     
      return
!-------------------end-of-vperm--------------------------------------- 
!-----------------------------------------------------------------------
      end

      subroutine dperm (nrow,a,ja,ia,ao,jao,iao,perm,qperm,job)
      integer nrow,ja(*),ia(nrow+1),jao(*),iao(nrow+1),perm(nrow),
     *qperm(*),job
      real*8 a(*),ao(*)
!-----------------------------------------------------------------------
! This routine permutes the rows and columns of a matrix stored in CSR
! format. i.e., it computes P A Q, where P, Q are permutation matrices.
! P maps row i into row perm(i) and Q maps column j into column qperm(j):
!      a(i,j)    becomes   a(perm(i),qperm(j)) in new matrix
! In the particular case where Q is the transpose of P (symmetric
! permutation of A) then qperm is not needed.
! note that qperm should be of length ncol (number of columns) but this
! is not checked.
!-----------------------------------------------------------------------
! Y. Saad, Sep. 21 1989 / recoded Jan. 28 1991.
!-----------------------------------------------------------------------
! on entry:
!----------
! n       = dimension of the matrix
! a, ja,
!    ia = input matrix in a, ja, ia format
! perm       = integer array of length n containing the permutation arrays
!        for the rows: perm(i) is the destination of row i in the
!         permuted matrix -- also the destination of column i in case
!         permutation is symmetric (job .le. 2)
!
! qperm      = same thing for the columns. This should be provided only
!         if job=3 or job=4, i.e., only in the case of a nonsymmetric
!        permutation of rows and columns. Otherwise qperm is a dummy
!
! job      = integer indicating the work to be done:
! * job = 1,2 permutation is symmetric  Ao :== P * A * transp(P)
!             job = 1      permute a, ja, ia into ao, jao, iao
!             job = 2 permute matrix ignoring real values.
! * job = 3,4 permutation is non-symmetric  Ao :== P * A * Q
!             job = 3      permute a, ja, ia into ao, jao, iao
!             job = 4 permute matrix ignoring real values.
!
! on return:
!-----------
! ao, jao, iao = input matrix in a, ja, ia format
!
! in case job .eq. 2 or job .eq. 4, a and ao are never referred to
! and can be dummy arguments.
! Notes:
!-------
!  1) algorithm is in place
!  2) column indices may not be sorted on return even  though they may be
!     on entry.
!----------------------------------------------------------------------c
! local variables
      integer locjob, mod
!
!     locjob indicates whether or not real values must be copied.
!
      locjob = mod(job,2)
!
! permute rows first
!
      call rperm (nrow,a,ja,ia,ao,jao,iao,perm,locjob)
!
! then permute columns
!
      locjob = 0
!
      if (job .le. 2) then
         call cperm (nrow,ao,jao,iao,ao,jao,iao,perm,locjob)
      else
         call cperm (nrow,ao,jao,iao,ao,jao,iao,qperm,locjob)
      endif
!
      return
!-------end-of-dperm----------------------------------------------------
      end


!-----------------------------------------------------------------------
      subroutine rperm (nrow,a,ja,ia,ao,jao,iao,perm,job)
      integer nrow,ja(*),ia(nrow+1),jao(*),iao(nrow+1),perm(nrow),job
      real*8 a(*),ao(*)
!-----------------------------------------------------------------------
! this subroutine permutes the rows of a matrix in CSR format.
! rperm  computes B = P A  where P is a permutation matrix.
! the permutation P is defined through the array perm: for each j,
! perm(j) represents the destination row number of row number j.
! Youcef Saad -- recoded Jan 28, 1991.
!-----------------------------------------------------------------------
! on entry:
!----------
! n       = dimension of the matrix
! a, ja, ia = input matrix in csr format
! perm       = integer array of length nrow containing the permutation arrays
!        for the rows: perm(i) is the destination of row i in the
!         permuted matrix.
!         ---> a(i,j) in the original matrix becomes a(perm(i),j)
!         in the output  matrix.
!
! job      = integer indicating the work to be done:
!             job = 1      permute a, ja, ia into ao, jao, iao
!                       (including the copying of real values ao and
!                       the array iao).
!             job .ne. 1 :  ignore real values.
!                     (in which case arrays a and ao are not needed nor
!                      used).
!
!------------
! on return:
!------------
! ao, jao, iao = input matrix in a, ja, ia format
! note :
!        if (job.ne.1)  then the arrays a and ao are not used.
!----------------------------------------------------------------------c
!           Y. Saad, May  2, 1990                                      c
!----------------------------------------------------------------------c
      logical values
      values = (job .eq. 1)
!
!     determine pointers for output matix.
!
      do 50 j=1,nrow
         i = perm(j)
         iao(i+1) = ia(j+1) - ia(j)
 50   continue
!
! get pointers from lengths
!
      iao(1) = 1
      do 51 j=1,nrow
         iao(j+1)=iao(j+1)+iao(j)
 51   continue
!
! copying
!
      do 100 ii=1,nrow
!
! old row = ii  -- new row = iperm(ii) -- ko = new pointer
!
         ko = iao(perm(ii))
         do 60 k=ia(ii), ia(ii+1)-1
            jao(ko) = ja(k)
            if (values) ao(ko) = a(k)
            ko = ko+1
 60      continue
 100  continue
!
      return
!---------end-of-rperm -------------------------------------------------
!-----------------------------------------------------------------------
      end

!-----------------------------------------------------------------------
      subroutine cperm (nrow,a,ja,ia,ao,jao,iao,perm,job)
      integer nrow,ja(*),ia(nrow+1),jao(*),iao(nrow+1),perm(*), job
      real*8 a(*), ao(*)
!-----------------------------------------------------------------------
! this subroutine permutes the columns of a matrix a, ja, ia.
! the result is written in the output matrix  ao, jao, iao.
! cperm computes B = A P, where  P is a permutation matrix
! that maps column j into column perm(j), i.e., on return
!      a(i,j) becomes a(i,perm(j)) in new matrix
! Y. Saad, May 2, 1990 / modified Jan. 28, 1991.
!-----------------------------------------------------------------------
! on entry:
!----------
! nrow       = row dimension of the matrix
!
! a, ja, ia = input matrix in csr format.
!
! perm      = integer array of length ncol (number of columns of A
!         containing the permutation array  the columns:
!         a(i,j) in the original matrix becomes a(i,perm(j))
!         in the output matrix.
!
! job      = integer indicating the work to be done:
!             job = 1      permute a, ja, ia into ao, jao, iao
!                       (including the copying of real values ao and
!                       the array iao).
!             job .ne. 1 :  ignore real values ao and ignore iao.
!
!------------
! on return:
!------------
! ao, jao, iao = input matrix in a, ja, ia format (array ao not needed)
!
! Notes:
!-------
! 1. if job=1 then ao, iao are not used.
! 2. This routine is in place: ja, jao can be the same.
! 3. If the matrix is initially sorted (by increasing column number)
!    then ao,jao,iao  may not be on return.
!
!----------------------------------------------------------------------c
! local parameters:
      integer k, i, nnz
!
      nnz = ia(nrow+1)-1
      do 100 k=1,nnz
         jao(k) = perm(ja(k))
 100  continue
!
!     done with ja array. return if no need to touch values.
!
      if (job .ne. 1) return
!
! else get new pointers -- and copy values too.
!
      do 1 i=1, nrow+1
         iao(i) = ia(i)
 1    continue
!
      do 2 k=1, nnz
         ao(k) = a(k)
 2    continue
!
      return
!---------end-of-cperm--------------------------------------------------
!-----------------------------------------------------------------------
      end

!------ MULTI COLORING -------------------------------------------------
      subroutine multic (n,ja,ia,ncol,kolrs,il,iord,maxcol,ierr) 
      integer n, ja(*),ia(n+1),kolrs(n),iord(n),il(maxcol+1),ierr 
!-----------------------------------------------------------------------
!     multicoloring ordering -- greedy algorithm -- 
!     determines the coloring permutation and sets up
!     corresponding data structures for it.
!-----------------------------------------------------------------------
! on entry
! -------- 
! n     = row and column dimention of matrix
! ja    = column indices of nonzero elements of matrix, stored rowwise.
! ia    = pointer to beginning of each row in ja.
! maxcol= maximum number of colors allowed -- the size of il is
!         maxcol+1 at least. Note: the number of colors does not
!         exceed the maximum degree of each node +1.
! iord  = en entry iord gives the order of traversal of the nodes
!         in the multicoloring algorithm. If there is no preference 
!         then set iord(j)=j for j=1,...,n
!
! on return
! --------- 
! ncol  = number of colours found 
! kolrs = integer array containing the color number assigned to each node 
! il    = integer array containing the pointers to the
!         beginning of each color set. In the permuted matrix
!         the rows /columns il(kol) to il(kol+1)-1 have the same color.
! iord  = permutation array corresponding to the multicolor ordering.
!         row number i will become row nbumber iord(i) in permuted 
!         matrix. (iord = destination permutation array).
! ierr  = integer. Error message. normal return ierr = 0. If ierr .eq.1
!         then the array il was overfilled. 
! 
!-----------------------------------------------------------------------
!     
      integer kol, i, j, k, maxcol, mycol 
!     
      ierr = 0
      do 1 j=1, n
         kolrs(j) = 0
 1    continue 
      do 11 j=1, maxcol
         il(j) = 0
 11   continue
!     
      ncol = 0
!     
!     scan all nodes 
!     
      do 4 ii=1, n
         i = iord(ii) 
!     
!     look at adjacent nodes to determine colors already assigned
!     
         mcol = 0
         do 2 k=ia(i), ia(i+1)-1
            j = ja(k)
            icol = kolrs(j)
            if (icol .ne. 0) then
               mcol = max(mcol,icol) 
!     
!     il used as temporary to record already assigned colors.
!     
               il(icol) = 1 
            endif
 2       continue
!     
!     taken colors determined. scan il until a slot opens up.
!     
         mycol = 1
 3       if (il(mycol) .eq. 1) then
            mycol = mycol+1 
            if (mycol .gt. maxcol) goto 99
            if (mycol .le. mcol) goto 3
         endif
!     
!     reset il to zero for next nodes
!     
         do 35 j=1, mcol
            il(j) = 0
 35      continue
!     
!     assign color and update number of colors so far
!     
         kolrs(i) = mycol
         ncol = max(ncol,mycol)
 4    continue
!     
!     every node has now been colored. Count nodes of each color
!     
      do 6 j=1, n
         kol = kolrs(j)+1
         il(kol) = il(kol)+1
 6    continue
!     
!     set pointers il
!     
      il(1) = 1
      do 7 j=1, ncol
         il(j+1) = il(j)+il(j+1)
 7    continue
!     
!     set iord
!     
      do 8 j=1, n
         kol = kolrs(j) 
         iord(j) = il(kol)
         il(kol) = il(kol)+1
 8    continue
!     
!     shift il back 
!     
      do 9 j=ncol,1,-1
         il(j+1) = il(j)
 9    continue
      il(1) = 1
!     
      return
 99   ierr = 1
      return
!----end-of-multic------------------------------------------------------
!-----------------------------------------------------------------------
      end

      subroutine filter(n,job,drptol,a,ja,ia,b,jb,ib,len,ierr)
      real*8 a(*),b(*),drptol
      integer ja(*),jb(*),ia(*),ib(*),n,job,len,ierr
!-----------------------------------------------------------------------
!     This module removes any elements whose absolute value
!     is small from an input matrix A and puts the resulting
!     matrix in B.  The input parameter job selects a definition
!     of small.
!-----------------------------------------------------------------------
! on entry:
!---------
!  n	 = integer. row dimension of matrix
!  job   = integer. used to determine strategy chosen by caller to
!         drop elements from matrix A. 
!          job = 1  
!              Elements whose absolute value is less than the
!              drop tolerance are removed.
!          job = 2
!              Elements whose absolute value is less than the 
!              product of the drop tolerance and the Euclidean
!              norm of the row are removed. 
!          job = 3
!              Elements whose absolute value is less that the
!              product of the drop tolerance and the largest
!              element in the row are removed.
! 
! drptol = real. drop tolerance used for dropping strategy.
! a	
! ja
! ia     = input matrix in compressed sparse format
! len	 = integer. the amount of space available in arrays b and jb.
!
! on return:
!---------- 
! b	
! jb
! ib    = resulting matrix in compressed sparse format. 
! 
! ierr	= integer. containing error message.
!         ierr .eq. 0 indicates normal return
!         ierr .gt. 0 indicates that there is'nt enough
!         space is a and ja to store the resulting matrix.
!         ierr then contains the row number where filter stopped.
! note:
!------ This module is in place. (b,jb,ib can ne the same as 
!       a, ja, ia in which case the result will be overwritten).
!----------------------------------------------------------------------c
!           contributed by David Day,  Sep 19, 1989.                   c
!----------------------------------------------------------------------c
! local variables
      real*8 norm,loctol
      integer index,row,k,k1,k2 
!
      index = 1
      do 10 row= 1,n
         k1 = ia(row)
         k2 = ia(row+1) - 1
         ib(row) = index
         goto (100,200,300) job
 100     norm = 1.0d0
         goto 400
 200     norm = 0.0d0
         do 22 k = k1,k2
            norm = norm + a(k) * a(k)
 22      continue
         norm = sqrt(norm)
         goto 400
 300     norm = 0.0d0
         do 23 k = k1,k2
            if( abs(a(k))  .gt. norm) then
               norm = abs(a(k))
            endif
 23      continue
 400     loctol = drptol * norm
         do 30 k = k1,k2
         if( abs(a(k)) .gt. loctol)then 
               if (index .gt. len) then
               ierr = row 
               return
            endif
            b(index) =  a(k)
            jb(index) = ja(k)
            index = index + 1
         endif
 30   continue
 10   continue
      ib(n+1) = index
      ierr = 0
      return
!--------------------end-of-filter -------------------------------------
!-----------------------------------------------------------------------
      end

!-----------------------------------------------------------------------
      subroutine getu (n,a,ja,ia,ao,jao,iao)
      integer n, ia(*), ja(*), iao(*), jao(*)
      real*8 a(*), ao(*)
!------------------------------------------------------------------------
! this subroutine extracts the upper triangular part of a matrix
! and writes the result ao, jao, iao. The routine is in place in
! that ao, jao, iao can be the same as a, ja, ia if desired.
!-----------
! on input:
!
! n     = dimension of the matrix a.
! a, ja,
!    ia = matrix stored in a, ja, ia, format
! On return:
! ao, jao,
!    iao = upper triangular matrix (upper part of a)
!      stored in compressed sparse row format
! note: the diagonal element is the last element in each row.
! i.e. in  a(ia(i+1)-1 )
! ao, jao, iao may be the same as a, ja, ia on entry -- in which case
! getu will overwrite the result on a, ja, ia.
!
!------------------------------------------------------------------------
! local variables
      real*8 t
      integer ko, k, i, kdiag, kfirst
      ko = 0
      do  7 i=1, n
         kfirst = ko+1
         kdiag = 0
         do 71 k = ia(i), ia(i+1) -1
            if (ja(k)  .lt. i) goto 71
            ko = ko+1
            ao(ko) = a(k)
            jao(ko) = ja(k)
            if (ja(k)  .eq. i) kdiag = ko
 71      continue
         if (kdiag .eq. 0 .or. kdiag .eq. kfirst) goto 72
!     exchange
         t = ao(kdiag)
         ao(kdiag) = ao(kfirst)
         ao(kfirst) = t
!
         k = jao(kdiag)
         jao(kdiag) = jao(kfirst)
         jao(kfirst) = k
 72      iao(i) = kfirst
 7    continue
!     redefine iao(n+1)
      iao(n+1) = ko+1
      return
!----------end-of-getu -------------------------------------------------
!-----------------------------------------------------------------------
      end

      subroutine submat(n,job,i1,i2,j1,j2,a,ja,ia,nr,nc,ao,jao,iao)
      integer n,job,i1,i2,j1,j2,nr,nc,ia(*),ja(*),jao(*),iao(*)
      real*8 a(*),ao(*)
!-----------------------------------------------------------------------
! extracts the submatrix A(i1:i2,j1:j2) and puts the result in
! matrix ao,iao,jao
!---- In place: ao,jao,iao may be the same as a,ja,ia.
!--------------
! on input
!---------
! n      = row dimension of the matrix
! i1,i2 = two integers with i2 .ge. i1 indicating the range of rows to be
!          extracted.
! j1,j2 = two integers with j2 .ge. j1 indicating the range of columns
!         to be extracted.
!         * There is no checking whether the input values for i1, i2, j1,
!           j2 are between 1 and n.
! a,
! ja,
! ia    = matrix in compressed sparse row format.
!
! job      = job indicator: if job .ne. 1 then the real values in a are NOT
!         extracted, only the column indices (i.e. data structure) are.
!         otherwise values as well as column indices are extracted...
!
! on output
!--------------
! nr      = number of rows of submatrix
! nc      = number of columns of submatrix
!        * if either of nr or nc is nonpositive the code will quit.
!
! ao,
! jao,iao = extracted matrix in general sparse format with jao containing
!      the column indices,and iao being the pointer to the beginning
!      of the row,in arrays a,ja.
!----------------------------------------------------------------------c
!           Y. Saad, Sep. 21 1989                                      c
!----------------------------------------------------------------------c
      nr = i2-i1+1
      nc = j2-j1+1
!
      if ( nr .le. 0 .or. nc .le. 0) return
!
      klen = 0
!
!     simple procedure. proceeds row-wise...
!
      do 100 i = 1,nr
         ii = i1+i-1
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         iao(i) = klen+1
!-----------------------------------------------------------------------
         do 60 k=k1,k2
            j = ja(k)
            if (j .ge. j1 .and. j .le. j2) then
               klen = klen+1
               if (job .eq. 1) ao(klen) = a(k)
               jao(klen) = j - j1+1
            endif
 60      continue
 100  continue
      iao(nr+1) = klen+1
      return
!------------end-of submat----------------------------------------------
!-----------------------------------------------------------------------
      end

!-----------------------------------------------------------------------
      subroutine infdia (n,ja,ia,ind,idiag) 
      integer ia(*), ind(*), ja(*)
!-----------------------------------------------------------------------
!     obtains information on the diagonals of A. 
!----------------------------------------------------------------------- 
! this subroutine finds the lengths of each of the 2*n-1 diagonals of A
! it also outputs the number of nonzero diagonals found. 
!----------------------------------------------------------------------- 
! on entry:
!---------- 
! n	= dimension of the matrix a.
!
! a,    ..... not needed here.
! ja, 			
! ia    = matrix stored in csr format
!
! on return:
!----------- 
!
! idiag = integer. number of nonzero diagonals found. 
! 
! ind   = integer array of length at least 2*n-1. The k-th entry in
!         ind contains the number of nonzero elements in the diagonal
!         number k, the numbering beeing from the lowermost diagonal
!         (bottom-left). In other words ind(k) = length of diagonal
!         whose offset wrt the main diagonal is = - n + k.
!----------------------------------------------------------------------c
!           Y. Saad, Sep. 21 1989                                      c
!----------------------------------------------------------------------c
      n2= n+n-1
      do 1 i=1,n2
         ind(i) = 0
 1    continue
      do 3 i=1, n
         do 2 k=ia(i),ia(i+1)-1
            j = ja(k)
            ind(n+j-i) = ind(n+j-i) +1
 2       continue 
 3    continue
!     count the nonzero ones.
      idiag = 0 
      do 41 k=1, n2
         if (ind(k) .ne. 0) idiag = idiag+1
 41   continue
      return
! done
!------end-of-infdia ---------------------------------------------------
!-----------------------------------------------------------------------
      end

      subroutine csrdia (n,idiag,job,a,ja,ia,ndiag,diag,ioff,ao,jao,iao,
     *ind)
      real*8 diag(ndiag,idiag), a(*), ao(*)
      integer ia(*), ind(*), ja(*), jao(*), iao(*), ioff(*)
!----------------------------------------------------------------------- 
! Compressed sparse row     to    diagonal format
!----------------------------------------------------------------------- 
! this subroutine extracts  idiag diagonals  from the  input matrix a,
! a, ia, and puts the rest of  the matrix  in the  output matrix ao,
! jao, iao.  The diagonals to be extracted depend  on the  value of job
! (see below for details.)  In  the first  case, the  diagonals to be
! extracted are simply identified by  their offsets  provided in ioff
! by the caller.  In the second case, the  code internally determines
! the idiag most significant diagonals, i.e., those  diagonals of the
! matrix which  have  the  largest  number  of  nonzero elements, and
! extracts them.
!----------------------------------------------------------------------- 
! on entry:
!---------- 
! n	= dimension of the matrix a.
! idiag = integer equal to the number of diagonals to be extracted. 
!         Note: on return idiag may be modified.
! a, ja, 			
!    ia = matrix stored in a, ja, ia, format
! job	= integer. serves as a job indicator.  Job is better thought 
!         of as a two-digit number job=xy. If the first (x) digit
!         is one on entry then the diagonals to be extracted are 
!         internally determined. In this case csrdia exctracts the
!         idiag most important diagonals, i.e. those having the largest
!         number on nonzero elements. If the first digit is zero
!         then csrdia assumes that ioff(*) contains the offsets 
!         of the diagonals to be extracted. there is no verification 
!         that ioff(*) contains valid entries.
!         The second (y) digit of job determines whether or not
!         the remainder of the matrix is to be written on ao,jao,iao.
!         If it is zero  then ao, jao, iao is not filled, i.e., 
!         the diagonals are found  and put in array diag and the rest is
!         is discarded. if it is one, ao, jao, iao contains matrix
!         of the remaining elements.
!         Thus:
!         job= 0 means do not select diagonals internally (pick those
!                defined by ioff) and do not fill ao,jao,iao
!         job= 1 means do not select diagonals internally 
!                      and fill ao,jao,iao
!         job=10 means  select diagonals internally 
!                      and do not fill ao,jao,iao
!         job=11 means select diagonals internally 
!                      and fill ao,jao,iao
! 
! ndiag = integer equal to the first dimension of array diag.
!
! on return:
!----------- 
!
! idiag = number of diagonals found. This may be smaller than its value 
!         on entry. 
! diag  = real array of size (ndiag x idiag) containing the diagonals
!         of A on return
!          
! ioff  = integer array of length idiag, containing the offsets of the
!   	  diagonals to be extracted.
! ao, jao
!  iao  = remainder of the matrix in a, ja, ia format.
! work arrays:
!------------ 
! ind   = integer array of length 2*n-1 used as integer work space.
!         needed only when job.ge.10 i.e., in case the diagonals are to
!         be selected internally.
!
! Notes:
!-------
!    1) The algorithm is in place: ao, jao, iao can be overwritten on 
!       a, ja, ia if desired 
!    2) When the code is required to select the diagonals (job .ge. 10) 
!       the selection of the diagonals is done from left to right 
!       as a result if several diagonals have the same weight (number 
!       of nonzero elemnts) the leftmost one is selected first.
!-----------------------------------------------------------------------
      job1 = job/10
      job2 = job-job1*10
      if (job1 .eq. 0) goto 50
      n2 = n+n-1
      call infdia(n,ja,ia,ind,idum)
!----------- determine diagonals to  accept.---------------------------- 
!----------------------------------------------------------------------- 
      ii = 0
 4    ii=ii+1
      jmax = 0
      do 41 k=1, n2
         j = ind(k)
         if (j .le. jmax) goto 41
         i = k
         jmax = j
 41   continue
      if (jmax .le. 0) then
         ii = ii-1
         goto 42
      endif
      ioff(ii) = i-n
      ind(i) = - jmax
      if (ii .lt.  idiag) goto 4
 42   idiag = ii
!---------------- initialize diago to zero ----------------------------- 
 50   continue
      do 55 j=1,idiag
         do 54 i=1,n
            diag(i,j) = 0.0d0
 54      continue
 55   continue
!----------------------------------------------------------------------- 
      ko = 1
!----------------------------------------------------------------------- 
! extract diagonals and accumulate remaining matrix.
!----------------------------------------------------------------------- 
      do 6 i=1, n
         do 51 k=ia(i),ia(i+1)-1 
            j = ja(k)
            do 52 l=1,idiag
               if (j-i .ne. ioff(l)) goto 52
               diag(i,l) = a(k)
               goto 51
 52         continue
!--------------- append element not in any diagonal to ao,jao,iao ----- 
            if (job2 .eq. 0) goto 51
            ao(ko) = a(k)
            jao(ko) = j
            ko = ko+1
 51      continue
         if (job2 .ne. 0 ) ind(i+1) = ko
 6    continue
      if (job2 .eq. 0) return
!     finish with iao
      iao(1) = 1
      do 7 i=2,n+1
         iao(i) = ind(i)
 7    continue
      return
!----------- end of dcsrdia ---------------------------------------------
!-----------------------------------------------------------------------
      end

