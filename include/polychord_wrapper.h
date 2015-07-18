#ifndef COSMO_PP_POLYCHORD_WRAPPER_H
#define COSMO_PP_POLYCHORD_WRAPPER_H

#ifdef POLYCHORD_IFORT
#define POLYRUN polycwraprun
#else
#define POLYRUN polycwraprun
#endif

namespace poly
{
    extern "C"
    {
        void POLYRUN(int ndims, int nderived, int nlive, int num_repeats, bool do_clustering, int ncluster, int feedback, bool calculate_post, int sigma_post, double thin_post, int *prior_types, double *prior_mins, double *prior_maxs, int *speeds, char *base_dir, char* file_root, bool read_resume, bool write_resume, int update_resume, bool write_live, double (*loglike)(double *theta, double *phi), double *logz, double *errorz, double *ndead, double *nlike, double *logzpluslogp, int num_grades, int *grade_dims, double *grade_fracs);
    }

    void run(int ndims, int nderived, int nlive, int num_repeats, bool do_clustering, int ncluster, int feedback, bool calculate_post, int sigma_post, double thin_post, int* prior_types, double *prior_mins, double *prior_maxs, int *speeds, const char *base_dir, const char *file_root, bool read_resume, bool write_resume, int update_resume, bool write_live, double (*loglike)(double *theta, double *phi), double *logz, double *errorz, double *ndead, double *nlike, double *logzpluslogp, int num_grades, int *grade_dims, double *grade_fracs)
    {
        char *bd = const_cast<char*>(base_dir);
        char *fr = const_cast<char*>(file_root);
        POLYRUN(ndims, nderived, nlive, num_repeats, do_clustering, ncluster, feedback, calculate_post, sigma_post, thin_post, prior_types, prior_mins, prior_maxs, speeds, bd, fr, read_resume, write_resume, update_resume, write_live, loglike, logz, errorz, ndead, nlike, logzpluslogp, num_grades, grade_dims, grade_fracs);
    }
} // namespace poly

#endif

