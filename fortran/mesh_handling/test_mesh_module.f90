program test_physical_group
  use mesh_module
  use matrix_writer
  use, intrinsic :: iso_c_binding
  implicit none
  
  type(PhysicalGroupData), allocatable :: physicalGroups(:)
  integer(c_size_t), allocatable :: nodeTagsAll(:)
  real(c_double), allocatable :: coordAll(:)
  real(c_double), allocatable :: parametricCoordAll(:)
  integer :: ierr, i
  
  ! Call the subroutine to get physical group nodes
  call get_physical_group_nodes("donut2d_mesh.msh", physicalGroups, ierr)
  
  if (ierr == 0) then
     ! Print physical group data for verification
     do i = 1, size(physicalGroups)
        print *, 'Physical Group ', i, &
                 ' (dim=', physicalGroups(i)%dim, &
                 ', tag=', physicalGroups(i)%tag, '): ', &
                 size(physicalGroups(i)%nodeTags), ' nodes'
     end do
     
     ! Clean up physical group data
     do i = 1, size(physicalGroups)
        if (allocated(physicalGroups(i)%nodeTags)) &
             deallocate(physicalGroups(i)%nodeTags)
        if (allocated(physicalGroups(i)%coord)) &
             deallocate(physicalGroups(i)%coord)
     end do
     deallocate(physicalGroups)
  else
     print *, 'Error in get_physical_group_nodes: ', ierr
  end if
 
  ! Call the subroutine to get all nodes
  call get_all_nodes("donut2d_mesh.msh", nodeTagsAll, coordAll, parametricCoordAll, ierr)
  print *, "Shape of coordAll: ", shape(coordAll)  
  if (ierr == 0) then
     ! Print all nodes data for verification
     print *, 'All Nodes: ', size(nodeTagsAll), ' nodes'
     print *, 'Node Tags (first 5): ', nodeTagsAll(1:min(5, size(nodeTagsAll)))
     print *, 'Coordinates (first 5): ', coordAll(1:min(5, size(coordAll)))
     if (size(parametricCoordAll) > 0) then
        print *, 'Parametric Coordinates (first 5): ', &
                 parametricCoordAll(1:min(5, size(parametricCoordAll)))
     else
        print *, 'No Parametric Coordinates Available'
     end if
     
  call write_matrix(reshape(coordAll,[3,size(coordAll)/3]), 'coordAll')

     ! Clean up all nodes data
     if (allocated(nodeTagsAll)) deallocate(nodeTagsAll)
     if (allocated(coordAll)) deallocate(coordAll)
     if (allocated(parametricCoordAll)) deallocate(parametricCoordAll)
  else
     print *, 'Error in get_all_nodes: ', ierr
  end if


  !call write_vector(coordAll, 'coordAll')
end program test_physical_group
