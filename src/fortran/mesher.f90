module mesher
implicit none
! ------------------------------------------------------------------------------
! diego domenzain @ CSM 2021
!
! this module is about building a mesh.
! It is based on my own Matlab code in
!
!      alles/graph-alg/mesher/
!
! * 'n_g2m_'        :: get the number of nodes in the graph.
! * 'g2m_m2g'       :: 
! * 'neigh_mesh_'   :: 
! * 'neigh_graph_'  :: 
! * 'neigh_type_'   :: 
! * 'nIJ'           :: 
! * 'IJ_'           :: 
! *
!
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
subroutine n_g2m_(nx,nz,a,n_g2m)
  ! it is assumed 'a' is a matrix with entries of 0 and 1,
  ! 0 : a point of no interest
  ! 1 : a point of interest
  ! ----------------------------------------------------------------------------
  ! get the number of nodes in the graph of the mesh 'a'
  ! ----------------------------------------------------------------------------
  integer, intent(in)     :: nx,nz,a(nz*nx)
  integer, intent(in out) :: n_g2m
  
  integer :: ia
  
  n_g2m = 0
  
  do ia = 1,nx*nz
    if (a(ia)==1) then
      n_g2m = n_g2m + 1
    endif
  enddo
  
end subroutine n_g2m_
! ------------------------------------------------------------------------------
subroutine g2m_m2g(nx,nz,n_g2m,a,graph2mesh,mesh2graph)
  ! it is assumed 'a' is a matrix with entries of 0 and 1,
  ! 0 : a point of no interest
  ! 1 : a point of interest
  ! ----------------------------------------------------------------------------
  ! make two dictionaries,
  ! graph2mesh : indexes are graph nodes, entries are mesh nodes
  ! mesh2graph : indexes are mesh nodes, entries are graph nodes
  ! ----------------------------------------------------------------------------
  integer, intent(in)     :: nx,nz,n_g2m,a(nz*nx)
  integer, intent(in out) :: graph2mesh(n_g2m), mesh2graph(nz*nx)
  
  integer :: ia
  
  i_g2m = 0
  
  do ia = 1,nx*nz
    if (a(ia)==1) then
      i_g2m = i_g2m + 1
      graph2mesh(i_g2m) = ia
      mesh2graph(ia) = i_g2m
    endif
  enddo
  
end subroutine g2m_m2g
! ------------------------------------------------------------------------------
subroutine neigh_mesh_(nx,nz,n_g2m,a,graph2mesh,neigh_mesh)
  ! it is assumed 'a' is a matrix with entries of 0 and 1,
  ! 0 : a point of no interest
  ! 1 : a point of interest
  ! ----------------------------------------------------------------------------
  ! neigh_mesh : row indexes are graph nodes.
  !              row entries are neighbors of that node, in the mesh. 
  ! ----------------------------------------------------------------------------
  integer, intent(in)     :: nx,nz,n_g2m,a(nz*nx),graph2mesh(n_g2m)
  integer, intent(in out) :: neigh_mesh(n_g2m,4)
  
  integer :: i_up,i_do,i_ri,i_le,i_g2m
  
  i_up = 0
  i_do = 0
  i_ri = 0
  i_le = 0
  
  do i_g2m = 1,n_g2m
    
    i_ri = graph2mesh(i_g2m) + nz
    i_up = graph2mesh(i_g2m) - 1
    i_le = graph2mesh(i_g2m) - nz
    i_do = graph2mesh(i_g2m) + 1
    
    ! --------------------------------------------------------------------------
    ! left edge
    if (graph2mesh(i_g2m).le.nz) then
      if (a(i_ri)==1) then
        neigh_mesh(i_g2m,1) = i_ri
      endif
      
      if (i_up.ge.1) then
        if (a(i_up)==1) then
          neigh_mesh(i_g2m,2) = i_up
        endif
      endif
      if (i_do.le.nz) then
        if (a(i_do)==1) then
          neigh_mesh(i_g2m,4) = i_do
        endif
      endif
    ! --------------------------------------------------------------------------
    ! right edge
    else if (graph2mesh(i_g2m).gt.nz*(nx-1)) then
      if (a(i_le)==1) then
        neigh_mesh(i_g2m,3) = i_le
      endif
      
      if (i_up.ge.nz*(nx-1)+1) then
        if (a(i_up)==1) then
          neigh_mesh(i_g2m,2) = i_up
        endif
      endif
      if (i_do.le.(nz*nx)) then
        if (a(i_do)==1) then
          neigh_mesh(i_g2m,4) = i_do
        endif
      endif
    ! --------------------------------------------------------------------------
    ! bottom edge
    else if (mod(graph2mesh(i_g2m),nz)==0) then
      if (a(i_up)==1) then
        neigh_mesh(i_g2m,2) = i_up
      endif
      
      if (i_ri.le.(nz*nx)) then
        if (a(i_ri)==1) then
          neigh_mesh(i_g2m,1) = i_ri
        endif
      endif
      if (i_le.ge.nz) then
        if (a(i_le)==1) then
          neigh_mesh(i_g2m,3) = i_le
        endif
      endif
    ! --------------------------------------------------------------------------
    ! top edge
    else if (mod(graph2mesh(i_g2m),nz)==1) then
      if (a(i_do)==1) then
        neigh_mesh(i_g2m,4) = i_do
      endif
      
      if (i_ri.le.(nz*(nx-1)+1)) then
        if (a(i_ri)==1) then
          neigh_mesh(i_g2m,1) = i_ri
        endif
      endif
      if (i_le.ge.1) then
        if (a(i_le)==1) then
          neigh_mesh(i_g2m,3) = i_le
        endif
      endif
    ! --------------------------------------------------------------------------
    ! inner nodes
    else
      if (a(i_ri)==1) then 
        neigh_mesh(i_g2m,1) = i_ri
      endif
      if (a(i_up)==1) then
        neigh_mesh(i_g2m,2) = i_up
      endif
      if (a(i_le)==1) then
        neigh_mesh(i_g2m,3) = i_le
      endif
      if (a(i_do)==1) then
        neigh_mesh(i_g2m,4) = i_do
      endif
    endif
  enddo
  
end subroutine neigh_mesh_
! ------------------------------------------------------------------------------
subroutine neigh_graph_(nx,nz,n_g2m,neigh_mesh,mesh2graph,neigh_graph)
  ! neigh_graph : row indexes are graph nodes.
  !               row entries are neighbors of that node, in the graph.
  ! ----------------------------------------------------------------------------
  integer, intent(in)     :: nx,nz,n_g2m,neigh_mesh(n_g2m,4),mesh2graph(nz*nx)
  integer, intent(in out) :: neigh_graph(n_g2m,4)
  
  integer :: i_g2m,i_nei
  
  do i_g2m = 1,n_g2m
    do i_nei = 1,4
      if (neigh_mesh(i_g2m,i_nei) ~= 0) then
        neigh_graph(i_g2m,i_nei) = mesh2graph(neigh_mesh(i_g2m,i_nei))
      endif
    enddo
  enddo
end subroutine neigh_graph_
! ------------------------------------------------------------------------------
subroutine  neigh_type_(nx,nz,n_g2m,a,graph2mesh,neigh_type)
  ! each node has a special type in a mesh-grid.
  ! neighbors that are (in the mesh):
  ! zero, non-zero, and next to the limits of the mesh.
  ! 
  ! In the case the mesh is a slice of the earth, and we are doing this for a PDE,
  ! this translates to:
  ! 
  ! neighbors that are non-zero = the pde (inner-node)
  ! neighbors that are zero = neumann bc
  ! neighbors that are next to the limits of the mesh = robin bc
  ! 
  ! neigh_type : row indexes are graph nodes.
  !              row entries are the type of neighbor for that node.
  ! 
  ! we define : (type,BC) = (1,inner) (-1,neumann) (0,robin)
  ! ----------------------------------------------------------------------------
  integer, intent(in)     :: nx,nz,n_g2m,a(nz*nx),graph2mesh(n_g2m)
  integer, intent(in out) :: neigh_type(n_g2m,4)
  
  integer :: inner,neuma,i_up,i_do,i_ri,i_le,i_g2m,i_nei
  
  do i_g2m = 1,n_g2m
    do i_nei = 1,4
      neigh_type(i_g2m,i_nei)=0
    enddo
  enddo

  inner =  1
  neuma = -1

  i_up = 0
  i_do = 0
  i_ri = 0
  i_le = 0

  do i_g2m = 1,n_g2m
    
    i_ri = graph2mesh(i_g2m) + nz
    i_up = graph2mesh(i_g2m) - 1
    i_le = graph2mesh(i_g2m) - nz
    i_do = graph2mesh(i_g2m) + 1
    ! --------------------------------------------------------------------------
    ! left edge
    if (graph2mesh(i_g2m).le.nz) then
      if (a(i_ri)==1) then
        neigh_type(i_g2m,1) = inner
      else
        neigh_type(i_g2m,1) = neuma
      endif
      
      if (i_up.ge.1) then
        if (a(i_up)==1) then
          neigh_type(i_g2m,2) = inner
        else
          neigh_type(i_g2m,2) = neuma
        endif
      endif
      if (i_do.le.nz) then
        if (a(i_do)==1) then
          neigh_type(i_g2m,4) = inner
        else
          neigh_type(i_g2m,4) = neuma
        endif
      endif
    ! --------------------------------------------------------------------------
    ! right edge
    elseif (graph2mesh(i_g2m).gt.nz*(nx-1))
      if (a(i_le)==1) then
        neigh_type(i_g2m,3) = inner
      else
        neigh_type(i_g2m,3) = neuma
      endif
      
      if (i_up.ge.nz*(nx-1)+1) then
        if (a(i_up)==1) then
          neigh_type(i_g2m,2) = inner
        else
          neigh_type(i_g2m,2) = neuma
        endif
      endif
      if (i_do.le.(nz*nx)) then
        if (a(i_do)==1) then
          neigh_type(i_g2m,4) = inner
        else
          neigh_type(i_g2m,4) = neuma
        endif
      endif
    ! --------------------------------------------------------------------------
    ! bottom edge
    elseif (mod(graph2mesh(i_g2m),nz)==0) then
      if (a(i_up)==1) then
        neigh_type(i_g2m,2) = inner
      else
        neigh_type(i_g2m,2) = neuma
      endif
      
      if (i_ri.le.(nz*nx)) then
        if (a(i_ri)==1) then
          neigh_type(i_g2m,1) = inner
        else
          neigh_type(i_g2m,1) = neuma
        endif
      endif
      if (i_le.ge.nz) then
        if (a(i_le)==1) then
          neigh_type(i_g2m,3) = inner
        else
          neigh_type(i_g2m,3) = neuma
        endif
      endif
    ! --------------------------------------------------------------------------
    ! top edge
    elseif (mod(graph2mesh(i_g2m),nz)==1) then
      if (a(i_do)==1) then
        neigh_type(i_g2m,4) = inner
      else
        neigh_type(i_g2m,4) = neuma
      endif
      
      if (i_ri.le.(nz*(nx-1)+1)) then
        if (a(i_ri)==1) then
          neigh_type(i_g2m,1) = inner
        else
          neigh_type(i_g2m,1) = neuma
        endif
      endif
      if (i_le.ge.1) then
        if (a(i_le)==1) then
          neigh_type(i_g2m,3) = inner
        else
          neigh_type(i_g2m,3) = neuma
        endif
      endif
    ! --------------------------------------------------------------------------
    ! inner nodes
    else
      if (a(i_ri)==1) then
        neigh_type(i_g2m,1) = inner
      else
        neigh_type(i_g2m,1) = neuma
      endif
      if (a(i_up)==1) then
        neigh_type(i_g2m,2) = inner
      else
        neigh_type(i_g2m,2) = neuma
      endif
      if (a(i_le)==1) then
        neigh_type(i_g2m,3) = inner
      else
        neigh_type(i_g2m,3) = neuma
      endif
      if (a(i_do)==1) then
        neigh_type(i_g2m,4) = inner
      else
        neigh_type(i_g2m,4) = neuma
      endif
    endif
  enddo
end subroutine neigh_type_
! ------------------------------------------------------------------------------
subroutine nIJ(n_g2m,neigh_type,n_ij,n_IJ)
  ! n_IJ : total number of non-zero entries of L (also length of I and J).
  ! n_ij : holds the info of how many entries in I(and J) belong to each node i.
  !        it is an array of size n_g2m by 1.
  ! ----------------------------------------------------------------------------
  integer, intent(in)     :: n_g2m,neigh_type(n_g2m,4)
  integer, intent(in out) :: n_ij(n_g2m),n_IJ

  integer :: i_g2m,i_nei

  do i_g2m=1,n_g2m
    n_i = 0
    do i_nei=1,4
      if (neigh_type(i_g2m,i_nei) .gt. 0) then
        n_i = n_i + neigh_type(i_g2m,i_nei)
      endif
    enddo
    n_ij(i_g2m) = n_i
    n_IJ = n_IJ + (n_i + 1)
  enddo
end subroutine nIJ
! ------------------------------------------------------------------------------
subroutine IJ_(n_g2m,n_ij,n_IJ,neigh_graph,I,J)
  integer, intent(in)     :: n_g2m,n_ij(n_g2m),n_IJ,neigh_graph(n_g2m,4)
  integer, intent(in out) :: I(n_IJ),J(n_IJ)
  
  integer :: i_g2m,ii,il,il_
  
  il = 1
  il_= 0
  
  do i_g2m=1,n_g2m
    ! -- end of line
    il_ = il_ + n_ij(i_g2m)+1
    ! -- in line
    ! I
    do ii=il,il_
      I(ii) = i_g2m
    enddo
    ! J
    J(il) = i_g2m
    i_nei = 1
    ii    = il+1
    do while (i_nei.le.4)
      if (neigh_graph(i_g2m,i_nei) .ne. 0) then
        J(ii) = neigh_graph(i_g2m,i_nei)
        ii = ii+1
      endif
      i_nei = i_nei + 1
    enddo
    ! -- begining of next line
    il = il_ + 1
  enddo
end subroutine IJ_
! ------------------------------------------------------------------------------
end module mesher