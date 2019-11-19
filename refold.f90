  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !
  !-----------------------------------------------------------------------
  subroutine refold ( ngm_g, mill_g, itoj, jtoi )
  !----------------------------------------------------------------------
  !
  !   Map the indices of G+G_0 into those of G 
  !   this is used to calculate electron-phonon matrix elements by
  !   refolding the k+q points into the first BZ (original k grid)
  !
  !   No parallelization on G-vecs at the moment  
  !   (actually this is done on the global array, but in elphel2.f90
  !   every processor has just a chunk of the array, I may need some
  !   communication)
  !
  !   No ultrasoft now
  !
  !   I use the rule : if not found then gmap = 0 
  !   Note that the map will be used only up to npwx (small sphere), 
  !   while the G-vectors lost in the process are on the surface of 
  !   the large sphere (density set).
  !
  !-----------------------------------------------------------------
  USE fft_base,      ONLY : dfftp
  USE io_global,     ONLY : stdout, meta_ionode
  USE io_epw,        ONLY : iukgmap
! SP: Sucidal. Produce too much data. Only use for debugging. 
!  USE control_flags, ONLY : iverbosity
  USE kfold
  !
  IMPLICIT NONE
  !
  INTEGER :: ngm_g
  !! Counter on G-vectors
  INTEGER :: mill_g(3,ngm_g)
  !!  Array of Miller indices of G-vectors in increasing order of G^2
  INTEGER :: jtoi(ngm_g)
  !! For the i-th G-vector in the sorted list, jtoi(i)
  !! returns its index in the unsorted list
  INTEGER :: itoj(ngm_g)
  !! itoj(i) returns the index of the G-vector in the sorted list
  !! that was at i-th position in the unsorted list
  INTEGER :: ig0 
  !! Counter on G_0 vectors
  INTEGER :: ig1
  !! Counter on G vectors
  INTEGER :: i, j, k
  !! Miller indices for G+G_0 vector
  INTEGER :: ig1_use, ig2_use
  !! Temporary G-vectors indices
  INTEGER :: indold, indnew
  !! Counter on G_0 vectors indices for writing to file
  !
  INTEGER :: ni, nj, nk
  INTEGER, ALLOCATABLE :: gvector2index(:, :, :)
  !! Mapping array of a G vector to its index in the unsorted list
  
  ALLOCATE( gmap(ngm_g,ng0vec) )
  gmap(:,:) = 0

  ! FHJ: Allocate mapping arrays to quickly search for G vectors.
  !
  ! Max miller indices (same convention as in module stick_set)
  ni = (dfftp%nr1-1)/2
  nj = (dfftp%nr2-1)/2
  nk = (dfftp%nr3-1)/2
  ALLOCATE( gvector2index(-ni:ni,-nj:nj,-nk:nk) )
  !
  ! FHJ: we don`t need to use itoj in the mapping arrays since mill_g already
  ! refers to the sorted G vectors.
  gvector2index(:,:,:) = 0
  DO ig2_use = 1, ngm_g
    i = mill_g(1,ig2_use)
    j = mill_g(2,ig2_use)
    k = mill_g(3,ig2_use)
    gvector2index(i,j,k) = ig2_use
  ENDDO

  !
  !  Loop on the inequivalent G_0 vectors
  !
  DO ig0 = 1, ng0vec
    !
    IF (ig0 .eq. 1) THEN
      WRITE(stdout,'(/5x,"Progress kgmap: ")',advance='no')
      indold = 0
    ENDIF
    indnew = nint( dble(ig0) / dble(ng0vec) * 40 )
    IF (indnew.ne.indold) WRITE(stdout,'(a)',advance='no') '#'
    indold = indnew
    !
    DO ig1 = 1, ngm_g
      !
      ig1_use = itoj(ig1)
      !
      !  the initial G vector
      !
      i = mill_g(1,ig1_use)
      j = mill_g(2,ig1_use)
      k = mill_g(3,ig1_use)
      !
      !  the final G+G_0 vector
      !
      i = i + g0vec_all(1,ig0)
      j = j + g0vec_all(2,ig0)
      k = k + g0vec_all(3,ig0)
      !
      ig2_use = 0
      ! FHJ: don`t change this for a cycle/goto, it can break vectorization!
      if (abs(i)<=ni .and. abs(j)<=nj .and. abs(k)<=nk) then
        ig2_use = gvector2index(i,j,k)
      endif
      gmap(ig1_use,ig0) = ig2_use
      !
    ENDDO
    !
  ENDDO
  !
  ! FHJ: deallocate G-vector mapping arrays
  DEALLOCATE( gvector2index )
  ! 
  !  output on file for electron-phonon matrix elements
  !
  IF (meta_ionode) then
    !FHJ: we now open *.kgmap as an unformatted file for performance
    DO ig1 = 1, ngm_g
      !WRITE(iukgmap,'(9i10)') (gmap(ig1,ig0), ig0 = 1, ng0vec)
      WRITE(iukgmap) gmap(ig1,1:ng0vec)
    ENDDO
    CLOSE(iukgmap)
  ENDIF
  WRITE(stdout,*)
  !
  RETURN
  !
  end subroutine refold

