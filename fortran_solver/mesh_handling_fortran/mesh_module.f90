module mesh_module
  use gmsh
  use, intrinsic :: iso_c_binding
  implicit none

  ! Define a derived type to store physical group data
  type :: PhysicalGroupData
     integer :: dim           ! Dimension of the physical group
     integer :: tag           ! Tag of the physical group
     integer(c_size_t), allocatable :: nodeTags(:)  ! Node tags
     real(c_double), allocatable :: coord(:)        ! Coordinates
  end type PhysicalGroupData

contains

  subroutine get_physical_group_nodes(filename, physicalGroups, ierr)
    character(len=*), intent(in) :: filename
    type(PhysicalGroupData), allocatable, intent(out) :: physicalGroups(:)
    integer, intent(out) :: ierr
    
    type(gmsh_t) :: g
    integer :: iPhysicalGroup, nPhysicalGroups
    integer(c_int), allocatable :: dimTags(:,:)
    
    ! Initialize GMSH
    call g%initialize(ierr=ierr)
    if (ierr /= 0) return
    
    ! Load the mesh file
    call g%open(filename, ierr)
    if (ierr /= 0) then
       call g%finalize(ierr)
       return
    end if
    
    ! Get all physical groups
    call g%model%getPhysicalGroups(dimTags, -1, ierr=ierr)
    if (ierr /= 0) then
       call g%finalize(ierr)
       return
    end if
    
    ! Allocate array to store data for each physical group
    nPhysicalGroups = size(dimTags, 2)
    allocate(physicalGroups(nPhysicalGroups))
    
    ! Loop through physical groups and store data
    do iPhysicalGroup = 1, nPhysicalGroups
        physicalGroups(iPhysicalGroup)%dim = dimTags(1, iPhysicalGroup)
        physicalGroups(iPhysicalGroup)%tag = dimTags(2, iPhysicalGroup)
        
        ! Get nodes for this physical group
        call g%model%mesh%getNodesForPhysicalGroup( &
             physicalGroups(iPhysicalGroup)%dim, &
             physicalGroups(iPhysicalGroup)%tag, &
             physicalGroups(iPhysicalGroup)%nodeTags, &
             physicalGroups(iPhysicalGroup)%coord, ierr)
        if (ierr /= 0) exit
    end do
    
    ! Finalize GMSH
    call g%finalize(ierr)
    
    ! Deallocate dimTags
    if (allocated(dimTags)) deallocate(dimTags)
  end subroutine get_physical_group_nodes

  subroutine get_nodes_connectivity(filename, connectivity, nodeTagsAll, coordAll, ierr)
    character(len=*), intent(in) :: filename
    integer(c_size_t), allocatable, intent(out) :: connectivity(:)
    integer(c_size_t), allocatable, intent(out) :: nodeTagsAll(:)
    real(c_double), allocatable, intent(out) :: coordAll(:)
    integer, intent(out) :: ierr
    
    real(c_double), allocatable :: parametricCoordAll(:)
    integer(c_size_t), allocatable :: elementTags(:)
    integer :: tag

    integer(c_int), allocatable :: elementTypes(:)    
    type(gmsh_t) :: g
    
    ! Initialize GMSH
    call g%initialize(ierr=ierr)
    if (ierr /= 0) return
    
    ! Load the mesh file
    call g%open(filename, ierr)
    if (ierr /= 0) then
       call g%finalize(ierr)
       return
    end if
   
    ! Get all element types
    call g%model%mesh%getElementTypes(elementTypes)

    ! Get connectivity
    tag = -1 ! Get all
    call g%model%mesh%getElementsByType(maxval(elementTypes), elementTags, connectivity, tag)
     
    ! Get all nodes and tags
    call g%model%mesh%getNodes(nodeTagsAll, coordAll, parametricCoordAll, ierr=ierr)
    
    ! Finalize GMSH
    call g%finalize(ierr)
  end subroutine get_nodes_connectivity

end module mesh_module
