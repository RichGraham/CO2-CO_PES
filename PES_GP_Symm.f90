module PES_details 
  interface PES 
     function PES(xStar) 
       implicit none 
       double precision:: PES
       double precision, dimension(:) ::  xStar
     end function PES
  end interface PES
end module PES_details

module GP_variables
  double precision, allocatable :: alpha (:), lScale(:), xTraining(:,:), xTrainingPerm(:,:)
  double precision expVar,NuggVar
  integer :: nDim=6
  integer :: nTraining=135 
end module GP_variables

  
program GP
  use GP_variables
  use PES_details
  implicit none
  
  double precision, allocatable::  xStar(:)
  double precision, allocatable::  gridData(:,:)
  integer :: nGrid=678 ! This is the raw size of the grid data
  integer i,nRMSE,j
  double precision :: dum, RMSE

  allocate (alpha(nTraining), lScale(nDim), xTraining(nDim,nTraining),xTrainingPerm(nDim,nTraining), xStar(nDim))
  allocate (gridData(nDim+1,nGrid))

  !====Load hyperparameters====
  open (unit = 7, file = "../HyperParams_Symm.dat")
  !Only need to read some as others are tied.
  read (7,*) lScale(4), lScale(3), lScale(2), lScale(1),expVar,NuggVar
  !Copy over the tied values
  lScale(5) = lScale(3)
  lScale(6) = lScale(4)
  print *,"HyperParams",lScale(1), lScale(2), lScale(3), lScale(4), lScale(5), lScale(6),expVar,NuggVar
  close(7)
  
  !====Load alpha coefficients====
  open (unit = 7, file = "../alpha_Symm.dat")
  do i=1,nTraining
     read (7,*) alpha(i)
     print *,"alpha ",i, alpha(i)
  end do
  close(7)


  !====Load training data x values ====
  open (unit = 7, file = "../xTraining.dat")
  do i=1,nTraining
     read (7,*) xTraining(1,i), xTraining(2,i), xTraining(3,i), xTraining(4,i), xTraining(5,i), xTraining(6,i)
     print *,"xTraining ", xTraining(1,i), xTraining(2,i), xTraining(3,i), xTraining(4,i), xTraining(5,i), xTraining(6,i)
  end do
  close(7)
  !! Permute the training vectors
  xTrainingPerm = xTraining
  do i=1,nTraining
     xTrainingPerm(3,i)=xTraining(5,i)
     xTrainingPerm(5,i)=xTraining(3,i)
     xTrainingPerm(4,i)=xTraining(6,i)
     xTrainingPerm(6,i)=xTraining(4,i)
  end do

  !====Load grid data ====
  open (unit = 7, file = "../rInv_HIt/co2CO_1074.lhc")
  !open (unit = 7, file = "../co2CO_6744.lhc")
  do i=1,nGrid
     read (7,*) gridData(1,i), gridData(2,i), gridData(3,i), gridData(4,i), gridData(5,i), gridData(6,i),gridData(7,i)
     !print *,"gridData ", gridData(1,i), gridData(2,i), gridData(3,i), gridData(4,i), gridData(5,i), gridData(6,i), gridData(7,i)
  end do
  close(7)


  
  xStar (1:6) = (/ 0.25054515,  0.30873283,  0.27150102,  0.31695355,  0.21821301,  0.26982914/)
  !xStar (1:6) = (/0.11845207,  0.11712603,  0.13018238,  0.13044285,  0.10768081,  0.10557117/)
    xStar (1:6) = (/0.47456526, 0.37298818, 0.49191578, 0.33957926, 0.3662692, 0.34482146 /)

  print *,"xStar", xStar
  dum=PES( xStar)
  print *,"Prediction =",PES( xStar)

  !====Compute RMSE====
  RMSE=0.0
  nRMSE = 0
  do i=1,nGrid
     if( gridData(7,i) < 0.005) then
        do j=1,nDim
           xStar(j) = gridData(j,i)
        end do
        RMSE  =  RMSE + (gridData(7,i)-PES(xStar))**2
        nRMSE = nRMSE +1
     end if
  end do
  print *,Sqrt(RMSE/(1.0*nRMSE))
  
end program GP
  
function PES(xStar)
  use GP_variables
  implicit none
  double precision, dimension(:) :: xStar
  double precision:: PES
  integer beta,i
  double precision kSqExp, kSqExpPerm, kKern

  kKern=0

  do i=1,nTraining
     kSqExp=1;
     kSqExpPerm=1;
     do beta=1,nDim
        kSqExp  =  kSqExp * ( exp( - (xStar(beta)-xTraining(beta,i))**2 /2.0/lScale(beta)**2) )
        kSqExpPerm  =  kSqExpPerm * ( exp( - (xStar(beta)-xTrainingPerm(beta,i))**2 /2.0/lScale(beta)**2) ) 
     end do
     kKern = kKern + alpha(i) * (kSqExp + kSqExpPerm)
  end do
  
  PES=kKern * expVar
end function PES