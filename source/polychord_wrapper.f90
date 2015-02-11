      module cpoly
          use iso_c_binding, only: c_funptr
          type(c_funptr) :: theloglike

          contains
              function loglike_f(theta, phi, context)
                  use iso_c_binding, only: c_double, c_f_procpointer
                  
                  implicit none
                  double precision, intent(in), dimension(:) :: theta
                  double precision, intent(out), dimension(:) :: phi
                  integer, intent(in) :: context
                  double precision :: loglike_f

                  interface
                      real(c_double) function loglike_proto(theta, phi, context)
                          use  iso_c_binding, only: c_int, c_double, c_ptr

                          implicit none

                          double precision theta, phi
                          !real(c_double), intent(in), dimension(:) :: theta
                          !real(c_double), intent(out), dimension(:) :: phi
                          integer(c_int), intent(in) :: context
                      end function loglike_proto
                  end interface

                  procedure(loglike_proto), pointer :: loglike_c

                  !print*,'fortran likelihood call'
                  !print*,'first arg = ',theta(1)
                  !print*,'second arg = ',theta(2)

                  call c_f_procpointer(theloglike,loglike_c)

                  loglike_f = loglike_c(theta(1), phi(1), context)
              end function loglike_f

              subroutine polycwraprun(ndims, nderived, nlive, num_repeats, &
                do_clustering, ncluster, feedback, calculate_post, &
                sigma_post, thin_post, prior_types, prior_mins, &
                prior_maxs, base_dir, file_root, &
                read_resume, write_resume, update_resume, write_live, &
                loglike, context, logz, errorz, ndead, nlike, &
                logzpluslogp) bind(c)
                      use iso_c_binding, only: c_int, c_bool, c_double, c_char, c_funptr, c_ptr, C_NULL_CHAR
                      use priors_module
                      use settings_module, only: program_settings, initialise_settings
                      use nested_sampling_module, only: NestedSampling

#ifdef COSMO_MPI
                      use mpi_module
#endif

                      implicit none

                      integer(c_int), intent(in), value :: ndims
                      integer(c_int), intent(in), value :: nderived
                      integer(c_int), intent(in), value :: nlive
                      integer(c_int), intent(in), value :: num_repeats
                      logical(c_bool), intent(in), value :: do_clustering
                      integer(c_int), intent(in), value :: ncluster
                      integer(c_int), intent(in), value :: feedback
                      logical(c_bool), intent(in), value :: calculate_post
                      integer(c_int), intent(in), value :: sigma_post
                      real(c_double), intent(in), value :: thin_post
                      integer(c_int), dimension(1), intent(in) :: prior_types
                      real(c_double), dimension(1), intent(in) :: prior_mins
                      real(c_double), dimension(1), intent(in) :: prior_maxs
                      character(kind=c_char,len=1), dimension(1), intent(in) :: base_dir
                      character(kind=c_char,len=1), dimension(1), intent(in) :: file_root
                      logical(c_bool), intent(in), value :: read_resume
                      logical(c_bool), intent(in), value :: write_resume
                      integer(c_int), intent(in), value :: update_resume
                      logical(c_bool), intent(in), value :: write_live
                      type(c_funptr), intent(in), value :: loglike
                      type(c_ptr), intent(in) :: context
                      double precision logz, errorz, ndead, nlike, logzpluslogp

                      character(len=100) :: froot, fdir
                      integer :: i

                      double precision, dimension(5) :: output_info
                      type(program_settings) :: settings
                      type(prior), allocatable, dimension(:) :: priors

                      double precision, dimension(1) :: minimums 
                      double precision, dimension(1) :: maximums
                      integer, dimension(1) :: hypercube_indices
                      integer, dimension(1) :: physical_indices

                      ! Temporary variables for initialising loglikelihoods
                      !double precision :: loglike1
                      !double precision, allocatable, dimension(:) :: theta
                      !double precision, allocatable, dimension(:) :: phi


                      froot = ' '
                      do i = 1, 100
                              if (file_root(i) == C_NULL_CHAR) then
                                      exit
                              else
                                      froot(i:i) = file_root(i)
                              end if
                      end do

                      fdir = ' '
                      do i = 1, 100
                              if (base_dir(i) == C_NULL_CHAR) then
                                      exit
                              else
                                      fdir(i:i) = base_dir(i)
                              end if
                      end do

                      theloglike = loglike

                      !print*,'PolyChord fortran run starting!'
                      !print*,'ndims = ',ndims
                      !print*,'nlinve = ',nlive
                      !print*,'num_repeats = ',num_repeats
                      !print*,'do_clustering = ',do_clustering
                      !print*,'ncluster = ',ncluster
                      !print*,'feedback = ',feedback
                      !print*,'calculate_post = ',calculate_post
                      !print*,'sigma_post = ',sigma_post
                      !print*,'thin_post = ',thin_post
                      !print*,'file_root = ',froot
                      !print*,'file_dir = ',fdir
                      !print*,'read_resume = ',read_resume
                      !print*,'write_resume = ',write_resume
                      !print*,'update_resume = ',update_resume
                      !print*,'write_live = ',write_live

                      !do i = 1, ndims
                      !print*,'parameter ',i,' has prior type ',prior_types(i)
                      !print*,'   min = ',prior_mins(i)
                      !print*,'   max = ',prior_maxs(i)
                      !end do

                      allocate(priors(ndims))
                      do i = 1, ndims
                                minimums(1) = prior_mins(i)
                                maximums(1) = prior_maxs(i)
                                hypercube_indices(1) = i
                                physical_indices(1) = i
                                if (prior_types(i) == 1) then
                                        call initialise_uniform(priors(i), &
                        hypercube_indices, physical_indices, minimums, maximums)
                                else if (prior_types(i) == 2) then
                                        call initialise_gaussian(priors(i), &
                        hypercube_indices, physical_indices, minimums, maximums)
                                else if (prior_types(i) == 3) then
                                        call initialise_log_uniform(priors(i), &
                        hypercube_indices, physical_indices, minimums, maximums)
                                else
                                        print*,'INVALID PRIOR TYPE ',i
                                        stop 1
                                end if
                      end do

                      settings%nDims = ndims
                      settings%nDerived = nderived
                      settings%nlive = nlive
                      settings%num_repeats = num_repeats
                      settings%do_clustering = logical(do_clustering)
                      settings%ncluster = ncluster
                      settings%feedback = feedback
                      settings%calculate_posterior = logical(calculate_post)
                      settings%sigma_posterior = sigma_post
                      settings%thin_posterior = thin_post
                      settings%file_root = froot
                      settings%base_dir = fdir
                      settings%read_resume = logical(read_resume)
                      settings%write_resume = logical(write_resume)
                      settings%update_resume = update_resume
                      settings%write_live = logical(write_live)

                      call initialise_settings(settings)   

                      ! ------- (1e) Initialise loglikelihood -----------------
                      ! This is only needed for a few things (e.g. generating a random correlated gaussian)
                      !allocate(theta(settings%nDims),phi(settings%nDerived))
                      !theta   = 0d0
                      !loglike1 = loglike_f(theta,phi,0)

#ifdef COSMO_MPI
                      output_info = NestedSampling(loglike_f, priors, &
                        settings, MPI_COMM_WORLD) 
#else
                      output_info = NestedSampling(loglike_f, priors, &
                        settings, 0) 
#endif

                      logz = output_info(1)
                      errorz = output_info(2)
                      ndead = output_info(3)
                      nlike = output_info(4)
                      logzpluslogp = output_info(5)

                      deallocate(priors)

              end subroutine
      end module
