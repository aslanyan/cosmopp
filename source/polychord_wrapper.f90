      module cpoly
          use iso_c_binding, only: c_funptr
          type(c_funptr) :: theloglike

          contains
              function loglike_f(theta, phi)
                  use iso_c_binding, only: c_double, c_f_procpointer
                  
                  implicit none
                  double precision, intent(in), dimension(:) :: theta
                  double precision, intent(out), dimension(:) :: phi
                  double precision :: loglike_f

                  interface
                      real(c_double) function loglike_proto(theta, phi)
                          use  iso_c_binding, only: c_int, c_double, c_ptr

                          implicit none

                          double precision theta, phi
                      end function loglike_proto
                  end interface

                  procedure(loglike_proto), pointer :: loglike_c

                  call c_f_procpointer(theloglike,loglike_c)

                  loglike_f = loglike_c(theta(1), phi(1))
              end function loglike_f

              subroutine polycwraprun(ndims, nderived, nlive, num_repeats, &
                do_clustering, ncluster, feedback, calculate_post, &
                sigma_post, thin_post, prior_types, prior_mins, &
                prior_maxs, base_dir, file_root, &
                read_resume, write_resume, update_resume, write_live, &
                loglike, logz, errorz, ndead, nlike, &
                logzpluslogp, num_grades, grade_dims, grade_fracs) bind(c)
                      use iso_c_binding, only: c_int, c_bool, c_double, c_char, c_funptr, c_ptr, C_NULL_CHAR
                      use ini_module, only:initialise_program
                      use params_module, only:add_parameter,param_type
                      use priors_module
                      use random_module, only:initialise_random
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
                      double precision logz, errorz, ndead, nlike, logzpluslogp
                      integer(c_int), intent(in), value :: num_grades
                      integer(c_int), dimension(1), intent(in) :: grade_dims
                      real(c_double), dimension(1), intent(in) :: grade_fracs

                      character(len=100) :: froot, fdir
                      integer :: i

                      double precision, dimension(5) :: output_info
                      type(program_settings) :: settings
                      type(prior), allocatable, dimension(:) :: priors
                      type(param_type),dimension(:),allocatable ::params
                      type(param_type),dimension(:),allocatable ::derived_params

                      double precision, dimension(1) :: minimums 
                      double precision, dimension(1) :: maximums
                      integer, dimension(1) :: hypercube_indices
                      integer, dimension(1) :: physical_indices

                      call initialise_random()

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

                      allocate(params(0),derived_params(0))

                      do i = 1, ndims
                                minimums(1) = prior_mins(i)
                                maximums(1) = prior_maxs(i)
                                hypercube_indices(1) = i
                                physical_indices(1) = i
                                if (prior_types(i) == 1) then
                                        call add_parameter(params,'p',&
                        'p',1,uniform_type,1,[prior_mins(i),prior_maxs(i)])
                                else if (prior_types(i) == 2) then
                                        call add_parameter(params,'p',&
                        'p',1,log_uniform_type,1,[prior_mins(i),prior_maxs(i)])
                                else if (prior_types(i) == 3) then
                                        call add_parameter(params,'p',&
                        'p',1,gaussian_type,1,[prior_mins(i),prior_maxs(i)])
                                else
                                        print*,'INVALID PRIOR TYPE ',i
                                        stop 1
                                end if
                      end do

                      settings%nlive = nlive
                      settings%num_repeats = num_repeats
                      settings%do_clustering = logical(do_clustering)
                      settings%base_dir = fdir
                      settings%file_root = froot
                      settings%write_resume = logical(write_resume)
                      settings%read_resume = logical(read_resume)
                      settings%write_live = logical(write_live)
                      settings%write_stats = .false.
                      settings%equals = .false.
                      settings%posteriors = logical(calculate_post)
                      settings%cluster_posteriors = .false.
                      settings%feedback = feedback
                      settings%update_resume = update_resume
                      settings%update_posterior = 1000
                      settings%boost_posterior = sigma_post
                      allocate(settings%grade_frac(num_grades))
                      allocate(settings%grade_dims(num_grades))

                      print *, num_grades
                      print *, size(settings%grade_dims)
                      do i = 1, num_grades
                            settings%grade_dims(i) = grade_dims(i)
                            settings%grade_frac(i) = grade_fracs(i)
                            print *, settings%grade_dims(i), ' ', settings%grade_frac(i)
                      end do

                      call initialise_program(settings,priors,params,derived_params)

                      print *, size(settings%grade_dims)

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
