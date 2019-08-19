  !LAST EDIT: Phil Maffettone 20
module likelihood
  use types
  use maths
  use class_proteins
  use ramachandran

  implicit none

  public

contains

  ! --------------------Public Contents-----------------------------------------
  ! FUNCTION: rama_likely(peptide_list, phi-psi-prby)
  ! SUBROUTINE: selection_hist(trip_prby,i_res,phi,psi,finalize)
  ! FUNCTION: rama_likely_wrtmax(peptide_list, phi-psi-prby)
  ! ----------------------------------------------------------------------------
  ! --------------------Private Contents----------------------------------------
  ! 
  ! ----------------------------------------------------------------------------

  function rama_likely(pep,trip_prby) result(r)
    ! FUNCTION: rama_likely(peptide_list, phi-psi-prby)
    ! Returns the log10 liklihood of the protein given the prior probability
    ! Depends on formality of Dunbrack ramachandran dist
    type(residue), intent(in) :: pep(:)                     !Protein
    real(dp), intent(in) :: trip_prby(:,:,:)                !Triplet prby as defined in ramachandran.f90
    real(dp) :: r                                           !Likelihood
    integer :: i,j,k                                        !Dummy integers

    r = 0.0_dp
    do i = 1, size(trip_prby,3)
       if (pep(i)%heteroatom) cycle
       j = floor((pep(i)%phi+180)/5)+1
       k = floor((pep(i)%psi+180)/5)+1
       j = 72*(j-1)+k
       r = r + log10(trip_prby(3,j,i))
    end do
  end function rama_likely
  ! -----------------------------------------------------------------------------

  function rama_likely_wrtmax(pep,trip_prby) result(r)
    ! FUNCTION: rama_likely_wrtmax(peptide_list, phi-psi-prby)
    ! Returns the log10 liklihood of the protein given the prior probability
    ! With respect to the maximum value of prior probability
    ! Depends on formality of Dunbrack ramachandran dist
    type(residue), intent(in) :: pep(:)                     !Protein
    real(dp), intent(in) :: trip_prby(:,:,:)                !Triplet prby as defined in ramachandran.f90
    real(dp) :: r                                           !Likelihood
    integer :: i,j,k                                        !Dummy integers

    r = 0.0_dp
    do i = 1, size(trip_prby,3)
       if (pep(i)%heteroatom) cycle
       j = floor((pep(i)%phi+180)/5)+1
       k = floor((pep(i)%psi+180)/5)+1
       j = 72*(j-1)+k
       r = r + log10(trip_prby(3,j,i))
    end do
    do i = 1, size(trip_prby,3)
       if (pep(i)%heteroatom) cycle
       r = r - log10(maxval(trip_prby(3,:,i)))
    end do
  end function rama_likely_wrtmax
  ! -----------------------------------------------------------------------------
  
  subroutine selection_hist(trip_prby,i_res,phi,psi,finalize)
    ! SUBROUTINE: selection_hist(phi-psi-prby,residue_selection,phi_sele,psi_sele,finalizer)
    real(dp),intent(inout) :: trip_prby(:,:,:)              !Triplet prby as defined in ramachandran.f90
    real(dp), intent(in) :: phi,psi                         !Angles selected
    real(dp) :: sum                                         !Sum for normaliztion
    integer, intent(in) :: i_res                            !Residue selected
    integer :: i                                            !Dummy index
    integer :: j                                            !index for phi/psi
    logical, intent(in), optional :: finalize               !Logical for final normalization


    if(present(finalize)) then
       if (finalize) then
          do i=1,size(trip_prby,3)
             sum = 0.
             do j=1,size(trip_prby,2)
                sum = sum+trip_prby(3,j,i)
             end do
             if (sum>0.) trip_prby(3,:,i) = trip_prby(3,:,i)/sum
          end do
          return
       end if
    end if

    j = 72*(floor(phi+180)/5) + (floor(psi+180)/5 + 1)
    trip_prby(3,j,i_res) = trip_prby(3,j,i_res) + 1.0_dp

  end subroutine selection_hist
  ! -----------------------------------------------------------------------------
end module likelihood
