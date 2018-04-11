module genericregression
  
  use iso_c_binding
  implicit none

  !
  ! Define an observation type that has a longitude and latitude coordinate
  ! with a value (d_i) and error (sigma_i).
  !
  type :: observation_t

     real(kind = c_double) :: value
     real(kind = c_double) :: sigma
     
     real(kind = c_double), dimension(1) :: lon
     real(kind = c_double), dimension(1) :: lat

  end type observation_t

  !
  ! Storage for the data to be loaded (dynamically allocatable)
  !
  integer :: nobs
  type(observation_t), dimension(:), allocatable :: obs

  real, parameter :: pi = 3.141592654
  
contains

  !
  ! This function is needed for interoperability between C++ and Fortran
  ! to convert and array of character to a string, see the open call in
  ! gvs2_loaddata_
  !
  function copystring(a) 
    character, intent(in) :: a(:)
    character(len = size(a)) :: copystring
    integer :: i
    do i = 1, size(a)
       copystring(i:i) = a(i)
    end do
  end function copystring

  !
  ! This function is called from the global state initialization constructor
  ! in the C++ code.
  !
  !  n         - The length of the character string filename
  !  filename  - The filename to load
  !  addobs    - A callback function to register points to evaluate for each
  !              observation
  !
  ! The addobs function must be called once and only once for each observation
  ! and includes all points needed to compute model predictions for
  ! the observation. For regression, this is a single point.
  !
  
  function gvs2_loaddata_(n, filename, addobs) bind(c)

    integer, intent(in) :: n
    character, dimension(n), intent(in) :: filename

    interface iaddobs
       function addobs(npoints, lons, lats) bind(c)
         use iso_c_binding
         integer :: addobs
         integer(kind = c_int), intent(in) :: npoints
         real(kind = c_double), dimension(npoints), intent(in) :: lons
         real(kind = c_double), dimension(npoints), intent(in) :: lats
       end function addobs
    end interface iaddobs
   
    integer :: gvs2_loaddata_
    
    !
    ! Local parameters
    !
    integer, parameter :: FUNIT = 10
    
    integer :: status
    integer :: i
    
    open(FUNIT, file = copystring(filename), status = 'old', iostat = status)
    if (status .ne. 0) then
       gvs2_loaddata_ = -1
       return
    end if
    
    !
    ! Custom text file format, first line contains the number of points
    ! Each subsequent line has floating point lon, lat, value, sigma
    !
    read(FUNIT, *) nobs
    
    !
    ! Allocate storage
    !
    allocate(obs(nobs))
    
    !
    ! Read each observation in
    !
    do i = 1, nobs
       
       read(FUNIT, *) obs(i)%lon(1), obs(i)%lat(1), obs(i)%value, obs(i)%sigma

       !
       ! Here we call addobs to register this observation with the inversion
       ! C++ code. For regression there is one point per observation and
       ! these are stored in the observation_t structure with length 1 arrays
       !
       if (addobs(1, obs(i)%lon, obs(i)%lat) .lt. 0) then
          gvs2_loaddata_ = -1
          return
       end if
       
    end do
    
    close(FUNIT)

    !
    ! Return 0 for success
    !
    gvs2_loaddata_ = 0
    
  end function gvs2_loaddata_

  !
  ! The compute prediction function is passed a single observation index and the values
  ! at the registered point(s) and needs to store the computed prediction in the
  ! prediction output parameter.
  !
  function gvs2_compute_prediction_(observationidx, npoints, values, unused, prediction) bind (c)
    use iso_c_binding
    
    integer, intent(in) :: observationidx
    integer, intent(in) :: npoints
    real(kind = c_double), dimension(npoints), intent(in) :: values
    real(kind = c_double), dimension(npoints), intent(out) :: unused
    real(kind = c_double), intent(out) :: prediction

    integer :: gvs2_compute_prediction_

    !
    ! For regression, there is no forward problem and the prediction is the value
    !
    prediction = values(1)

    gvs2_compute_prediction_ = 0

  end function gvs2_compute_prediction_

  !
  ! The compute likelihood is given all predictions and needs to compute the residuals,
  ! ie error between predictions and observations, and the negative log likelihood plus
  ! normalization term
  !
  function gvs2_compute_likelihood_(nobservation, hvalue, predictions, residuals, unused, like, norm) bind (c)
    use iso_c_binding
    
    integer, intent(in) :: nobservation
    real(kind = c_double), intent(in) :: hvalue
    real(kind = c_double), dimension(nobservation), intent(in) :: predictions
    real(kind = c_double), dimension(nobservation), intent(out) :: residuals
    real(kind = c_double), dimension(nobservation), intent(out) :: unused
    real(kind = c_double), intent(out) :: like
    real(kind = c_double), intent(out) :: norm

    integer :: gvs2_compute_likelihood_

    !
    ! Local parameters
    !
    real(kind = c_double) :: res
    integer :: i
    real(kind = c_double) :: n
    
    like = 0.0

    !
    ! First part of -ve log of norm, ie
    !       |
    !       v
    ! -log(1/sqrt(2 pi) * 1/(prod sigma_i))
    !
    norm = 0.5*log(2.0*PI)

    do i = 1, nobservation

       ! Loop through predictions and compute G(m)_i - d_i
       res = predictions(i) - obs(i)%value

       ! Store residual in output array
       residuals(i) = res

       ! Compute scaled data error using hierarchical term
       ! (in non-hierarchical, hvalue will be 1 and have no effect)
       n = hvalue * obs(i)%sigma

       !
       ! Add -ve log likelihood term (G(m)_i - d_i)^2/(2 sigma^2)
       !
       like = like + (res*res)/(2.0*n*n)

       !
       ! Add -ve log normalization term -log(1/sigma_i)
       !
       norm = norm + log(n)

    end do

    !
    ! Return 0 for success
    !
    gvs2_compute_likelihood_ = 0
    
  end function gvs2_compute_likelihood_

  !
  ! The save data function is used for synthetic observation generation
  ! The requirements for this function are to write to the specified file
  ! the observations previously loaded by the gvs2_loaddata_ function above,
  ! but with value overwritten by the values in the predictions array, and
  ! noiselevel outputted as the data errors.
  !
  function gvs2_savedata_(n, filename, noiselevel, nobservation, predictions) bind(c)

    integer, intent(in) :: n
    character, dimension(n), intent(in) :: filename

    real(kind = c_double), intent(in) :: noiselevel
    integer, intent(in) :: nobservation
    real(kind = c_double), dimension(nobservation), intent(in) :: predictions

    integer :: gvs2_savedata_

    !
    ! Local parameters
    !

    integer, parameter :: FUNIT = 10
    
    integer :: status
    integer :: i
    
    open(FUNIT, file = copystring(filename), status = 'replace', action = 'write', iostat = status)
    if (status .ne. 0) then
       gvs2_savedata_ = -1
       return
    end if
    
    write(FUNIT, *) nobservation
    
    do i = 1, nobservation
       
       write(FUNIT, *) obs(i)%lon(1), obs(i)%lat(1), predictions(i), noiselevel
       
    end do
    
    close(FUNIT)
    
    gvs2_savedata_ = 0
    
  end function gvs2_savedata_
  
end module genericregression
