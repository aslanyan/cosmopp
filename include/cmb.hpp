#ifndef COSMO_PP_CMB_HPP
#define COSMO_PP_CMB_HPP

#include <vector>

#include <table_function.hpp>
#include <cosmological_params.hpp>

struct precision;
struct background;
struct thermo;
struct perturbs;
struct bessels;
struct transfers;
struct primordial;
struct spectra;
struct nonlinear;
struct lensing;
struct output;

/// This class can be used to calculate cmb power spectra, transfer functions, etc. It uses the publicly availble code CLASS.
/// To do calculations one needs to first pre-initialize it, then initialize it, after which one can get the power spectra and the transfer functions.
/// Initialization can be done multiple times after pre-initialization with new parameters. The initialization is divided into two steps for faster calculations with different cosmological parameters.
/// Pre-initialization usually needs to be done once, for example in MCMC runs, then one can initialize multiple times leaving the pre-initialization parameters unchanged.
class CMB
{
public:
    /// Constructor.
    CMB() : preInit_(false), init_(false), primordialInitialize_(true), lensing_(false), includeTensors_(false) { allocate(); }

    /// Destructor.
    virtual ~CMB() { if(preInit_) preClean(); deAllocate(); }

    /// Pre-initialization routine. Always must be called after the constructor and before initialization can be done.
    /// \param lMax The maximum value of l for calculations. NOTE: Lensed cl-s will be calculated up to a value less than lMax, if requested. It is recommended to give lMax higher than requested lensed cl-s by 1000.
    /// \param wantAllL Defines if cl-s should be explicitly calculated for all l or not. By default the calculation is done for some l values, interpolated for the rest, for speedup. NOTE: If you're planning to use transfer functions make sure to turn on this flag.
    /// \param primordialInitialize If the primordial power spectrum is determined by only A_s and n_s this flag can be set to false. If the primordial power spectrum has some shape that requires calls to the function for each k, set this to true (true by default).
    /// \param includeTensors Specifies if tensor modes should be included in calculations.
    /// \param lMaxTensors Maximum l value for tensor modes. Only matters if includeTensors is true.
    /// \param kPerDecade Used only if primordialInitialize is true. Determines how many k points per decade need to be used for initializing the primordial power spectrum.
    /// \param kMin Used only if primordialInitialize is true. This is the lower end for initializing the primordial power spectrum.
    /// \param kMax Used only if primordialInitialize is true. This is the upper end for initializing the primordial power spectrum.
    virtual void preInitialize(int lMax, bool wantAllL = false, bool primordialInitialize = true, bool includeTensors = false, int lMaxTensors = 0, double kPerDecade = 100, double kMin = 1e-6, double kMax = 1.0);

    /// Initialization routine. Must be called after pre-initialization. Can be called multiple times in a row.
    /// \param params The cosmological parameters to use.
    /// \param wantT A flag specifying if the T mode (temperature) should be calculated.
    /// \param wantPol A flag specifying if polarization modes (E and B) should be calculated.
    /// \param wantLensing A flag specifying if lensing potential and lensed Cl-s should be calculated.
    /// \param wantMatterPs A flag specifying if the matter power spectrum should be calculated.
    /// \param zMaxPk The redshift up to which the matter power spectrum is needed. This parameter matters only if wantMatterPs is true.
    virtual void initialize(const CosmologicalParams& params, bool wantT = true, bool wantPol = false, bool wantLensing = false, bool wantMatterPs = false, double zMaxPk = 0);

    /// Retrieves the values of the calculated CMB power spectra.
    /// \param clTT A pointer to a vector where TT power spectra should be stored. Give NULL if not wanted. NOTE: Can only be requested if T modes have been calculated during initialization.
    /// \param clEE A pointer to a vector where EE power spectra should be stored. Give NULL if not wanted. NOTE: Can only be requested if polarization modes have been calculated during initialization.
    /// \param clTE A pointer to a vector where TE power spectra should be stored. Give NULL if not wanted. NOTE: Can only be requested if T modes and polarization modes have been calculated during initialization.
    /// \param clPP A pointer to a vector where PP power spectra (lensing potential) should be stored. Give NULL if not wanted. NOTE: Can only be requested if lensing has been calculated during initialization.
    /// \param clTP A pointer to a vector where TP power spectra should be stored. Give NULL if not wanted. NOTE: Can only be requested if T modes and lensing have been calculated during initialization.
    /// \param clEP A pointer to a vector where EP power spectra should be stored. Give NULL if not wanted. NOTE: Can only be requested if polarization modes and lensing have been calculated during initialization.
    /// \param clBB A pointer to a vector where BB power spectra should be stored. Give NULL if not wanted. NOTE: Can only be requested if polarization modes have been calculated during initialization.
    virtual void getCl(std::vector<double>* clTT, std::vector<double>* clEE = NULL, std::vector<double>* clTE = NULL, std::vector<double>* clPP = NULL, std::vector<double>* clTP = NULL, std::vector<double>* clEP = NULL, std::vector<double>* clBB = NULL);

    /// Retrieves the values of lensed Cl-s. This function can be called only if lensing has been calculated during initialization.
    /// \param clTT A pointer to a vector where lensed TT power spectra should be stored. Give NULL if not wanted. NOTE: Can only be requested if T modes have been calculated during initialization.
    /// \param clEE A pointer to a vector where lensed EE power spectra should be stored. Give NULL if not wanted. NOTE: Can only be requested if polarization modes have been calculated during initialization.
    /// \param clTE A pointer to a vector where lensed TE power spectra should be stored. Give NULL if not wanted. NOTE: Can only be requested if T and polarization modes have been calculated during initialization.
    /// \param clBB A pointer to a vector where lensed BB power spectra should be stored. Give NULL if not wanted. NOTE: Can only be requested if polarization modes have been calculated during initialization.
    virtual void getLensedCl(std::vector<double>* clTT, std::vector<double>* clEE = NULL, std::vector<double>* clTE = NULL, std::vector<double>* clBB = NULL);

    /// Retrieves CMB transfer functions. In order for this function to work properly all l calculation must be requested during pre-initialization.
    /// \param l The value of l for the transfer function.
    /// \param t A pointer to a map where the temperature transfer function should be stored. Give NULL if not wanted. NOTE: Can only be requested if T modes have been calculated during initialization.
    /// \param e A pointer to a map where the E mode transfer function should be stored. Give NULL if not wanted. NOTE: Can only be requested if polarization modes have been calculated during initialization.
    /// \param p A pointer to a map where the the lensing potential transfer function should be stored. Give NULL if not wanted. NOTE: Can only be requested if lensing has been calculated during initialization.
    virtual void getTransfer(int l, Math::TableFunction<double, double>* t, Math::TableFunction<double, double>* e = NULL, Math::TableFunction<double, double>* p = NULL);

    /// Retrieves the matter power spectrum. Should not be called unless wantMatterPs is set to true in initialize.
    /// \param z The redshift at which the matter power spectrum is wanted. Should not exceed zMaxPk of initialize.
    /// \param ps A pointer to a map where the matter power spectrum will be written.
    virtual void getMatterPs(double z, Math::TableFunction<double, double>* ps);

    /// Retrieves the matter transfer function. Should not be called unless wantMatterPs is set to true in initialize.
    /// \param z The redshift at which the matter transfer function is wanted. Should not exceed zMaxPk of initialize.
    /// \param tk A pointer to a map where the matter transfer function will be written.
    virtual void getMatterTransfer(double z, Math::TableFunction<double, double>* tk);

    /// Calculates sigma_8. Should not be called unless wantMatterPs is set to true in initialize.
    /// \return The sigma_8 value.
    virtual double sigma8();

    virtual void getLList(std::vector<double>& list) const;
    virtual void getLListLens(std::vector<double>& list) const;

protected:
    void preClean();
    void clean();

    void allocate();
    void deAllocate();

protected:
    const CosmologicalParams* params_;

    precision* pr_;
    background* br_;
    thermo* th_;
    perturbs* pt_;
    //bessels* bs_;
    transfers* tr_;
    primordial* pm_;
    spectra* sp_;
    nonlinear* nl_;
    lensing* le_;
    output* op_;

    bool primordialInitialize_;

    bool preInit_;
    bool init_;
    bool lensing_;
    int lMax_;
    
    bool includeTensors_;
    int lMaxTensors_;
    double kPerDecade_;
    double kMin_;
    double kMax_;

    double* dum1_[1000], *dum2_[1000];
};

#endif

