module genericregression
  
  use iso_c_binding
  implicit none
  
  type :: observation_t

     real(kind = c_double) :: ttime
     real(kind = c_double) :: sigma
     integer :: n
     
     real(kind = c_double), dimension(:), allocatable :: lon
     real(kind = c_double), dimension(:), allocatable :: lat
     real(kind = c_double), dimension(:), allocatable :: r
     real(kind = c_double), dimension(:), allocatable :: weight
     
  end type observation_t

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
  ! Compute the great circle distance between two lon, lat points at radius r
  ! using Haversine formula (accurate except for near antipodal points).
  !
  function greatcircledistkm(lon1, lat1, lon2, lat2, r)
    real(kind = c_double), intent(in) :: lon1
    real(kind = c_double), intent(in) :: lat1
    real(kind = c_double), intent(in) :: lon2
    real(kind = c_double), intent(in) :: lat2
    real(kind = c_double), intent(in) :: r

    real(kind = c_double) :: greatcircledistkm

    real(kind = c_double) :: phi1, theta1
    real(kind = c_double) :: phi2, theta2
    real(kind = c_double) :: hsinphi, hsintheta

    phi1 = (90.0 - lat1) * pi/180.0
    theta1 = lon1 * pi/180.0
    phi2 = (90.0 - lat2) * pi/180.0
    theta2 = lon2 * pi/180.0

    hsinphi = sin((phi1 - phi2)/2.0)
    hsintheta = sin((theta1 - theta2)/2.0)

    greatcircledistkm = 2.0*asin(sqrt(hsinphi*hsinphi + sin(phi1)*sin(phi2)*(hsintheta*hsintheta))) * r
    
  end function greatcircledistkm

  !
  ! Compute the cartesian distance between two lon, lat, r points
  ! by converting each to cartesian x, y, z points and computing
  ! the distance
  !
  function cartesiandistkm(lon1, lat1, r1, lon2, lat2, r2)
    real(kind = c_double), intent(in) :: lon1
    real(kind = c_double), intent(in) :: lat1
    real(kind = c_double), intent(in) :: r1
    real(kind = c_double), intent(in) :: lon2
    real(kind = c_double), intent(in) :: lat2
    real(kind = c_double), intent(in) :: r2

    real(kind = c_double) :: cartesiandistkm

    real(kind = c_double) :: phi1, theta1
    real(kind = c_double) :: phi2, theta2
    real(kind = c_double) :: x1, y1, z1
    real(kind = c_double) :: x2, y2, z2
    
    phi1 = (90.0 - lat1) * pi/180.0
    theta1 = lon1 * pi/180.0
    phi2 = (90.0 - lat2) * pi/180.0
    theta2 = lon2 * pi/180.0

    x1 = r1 * cos(theta1) * sin(phi1)
    y1 = r1 * sin(theta1) * sin(phi1)
    z1 = r1 * cos(phi1)
 
    x2 = r2 * cos(theta2) * sin(phi2)
    y2 = r2 * sin(theta2) * sin(phi2)
    z2 = r2 * cos(phi2)

    cartesiandistkm = sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)
    
  end function cartesiandistkm

  !
  ! Given a set of lon, lat, r points representing a path
  ! on or through a sphere, compute a distance weight for
  ! each point for computing travel time
  !
  subroutine computeweights(n, lons, lats, r, weights)

    integer, intent(in) :: n
    real(kind = c_double), dimension(n), intent(in) :: lons
    real(kind = c_double), dimension(n), intent(in) :: lats
    real(kind = c_double), dimension(n), intent(in) :: r
    real(kind = c_double), dimension(n), intent(out) :: weights

    integer :: i
    real(kind = c_double) :: dist
    
    weights = 0.0
    
    do i = 2, n

       if (r(i - 1) .eq. r(i)) then
          ! Great circle
          dist = greatcircledistkm(lons(i - 1), lats(i - 1), lons(i), lats(i), r(i))

          weights(i - 1) = weights(i - 1) + dist/2.0
          weights(i) = weights(i) + dist/2.0

       else
          ! Cartesian distance
          dist = cartesiandistkm(lons(i - 1), lats(i - 1), r(i - 1), lons(i), lats(i), r(i))

          weights(i - 1) = weights(i - 1) + dist/2.0
          weights(i) = weights(i) + dist/2.0
          
       end if

    end do

  end subroutine computeweights
          
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
    integer :: j
    
    open(FUNIT, file = copystring(filename), status = 'old', iostat = status)
    if (status .ne. 0) then
       gvs2_loaddata_ = -1
       return
    end if
    
    !
    ! First line contains the number of points
    !
    read(FUNIT, *) nobs
    
    !
    ! Allocate storage
    !
    allocate(obs(nobs))
    
    !
    ! Read each observation in. Each observation starts with a single line
    ! containing the observed travel time, the error and the number
    ! of path points. Then there is a separate line for each point
    ! consisting of lon, lat, radius.
    !
    do i = 1, nobs
       
       read(FUNIT, *) obs(i)%ttime, obs(i)%sigma, obs(i)%n

       allocate(obs(i)%lon(obs(i)%n))
       allocate(obs(i)%lat(obs(i)%n))
       allocate(obs(i)%r(obs(i)%n))
       allocate(obs(i)%weight(obs(i)%n))

       do j = 1, obs(i)%n

          read(FUNIT, *) obs(i)%lon(j), obs(i)%lat(j), obs(i)%r(j)

       end do

       !
       ! This subroutine computes the distance weights for each point
       !
       call computeweights(obs(i)%n, obs(i)%lon, obs(i)%lat, obs(i)%r, obs(i)%weight)

       !
       ! Here we call addobs to register this observation with the inversion
       ! C++ code. We pass in the lon and lat arrays we've just created.
       !
       if (addobs(obs(i)%n, obs(i)%lon, obs(i)%lat) .lt. 0) then
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
    ! Local
    !
    integer :: oi
    integer :: i
    real(kind = c_double) :: sum

    !
    ! For basic tomography, weights contains the distances and the model is velocity
    ! and is in the values input array so the predicted travel time is a sum over
    ! distance/velocity
    !

    !
    ! The passed observation index uses C index so is from 0 .. n - 1. Convert
    ! to fortran 1 .. n here.
    !
    oi = observationidx + 1
    sum = 0.0
    do i = 1, npoints

       sum = sum + obs(oi)%weight(i)/values(i)

    end do
       
    prediction = sum

    gvs2_compute_prediction_ = 0

  end function gvs2_compute_prediction_

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
       res = predictions(i) - obs(i)%ttime

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
   integer :: j

   open(FUNIT, file = copystring(filename), status = 'replace', action = 'write', iostat = status)
   if (status .ne. 0) then
      gvs2_savedata_ = -1
      return
   end if

   write(FUNIT, *) nobservation

   do i = 1, nobservation

      write(FUNIT, *) predictions(i), noiselevel, obs(i)%n

      do j = 1, obs(i)%n

         write(FUNIT, *) obs(i)%lon(j), obs(i)%lat(j), obs(i)%r(j)

      end do

   end do

   close(FUNIT)

   gvs2_savedata_ = 0
    
  end function gvs2_savedata_
    
end module genericregression
