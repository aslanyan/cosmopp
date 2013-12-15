#ifndef COSMO_PP_MASTER_HPP
#define COSMO_PP_MASTER_HPP

#include <vector>
#include <map>
#include <string>

#include <macros.hpp>

#include <gmd.h>

#include <healpix_map.h>

/// The MASTER algorithm for calculating the power spectrum.

/// This class makes use of the MASTER algorithm (astro-ph/0105302) for calculating the power spectrum.
class Master
{
public:
    /// Constructor.
    /// \param mask A reference to the mask (can be weights too).
    /// \param couplingKernelFileName The name of the file containing the coupling kernel of the mask. Can be calculated using calculateCouplingKernel.
    /// \param beam The beam function (including the pixel window function). The index is l.
    /// \param bins Specify this if you want to calculate the power spectrum on bins. Give NULL otherwise and it will calculate it for all l. The vector should contain the l values of the starting points of the bins. The last value will be an upper limit for the last bin. If this parameter is not NULL, the parameter lMax needs to be 0.
    /// \param lMax The maximum value of l up to which the power spectrum should be calculated. If bins are specified this value should be left to its default value of 0.
    Master(const Healpix_Map<double>& mask, const char* couplingKernelFileName, const std::vector<double>& beam, const std::vector<int>* bins, int lMax = 0);

    /// Constructor.
    /// \param maskName The name of the fits file containing the mask.
    /// \param couplingKernelFileName The name of the file containing the coupling kernel of the mask. Can be calculated using calculateCouplingKernel.
    /// \param beam The beam function (including the pixel window function). The index is l.
    /// \param bins Specify this if you want to calculate the power spectrum on bins. Give NULL otherwise and it will calculate it for all l. The vector should contain the l values of the starting points of the bins. The last value will be an upper limit for the last bin. If this parameter is not NULL, the parameter lMax needs to be 0.
    /// \param lMax The maximum value of l up to which the power spectrum should be calculated. If bins are specified this value should be left to its default value of 0.
    Master(const char* maskName, const char* couplingKernelFileName, const std::vector<double>& beam, const std::vector<int>* bins, int lMax = 0);

    /// Calculate the power spectrum for a given map.
    /// \param map The map to be used.
    void calculate(const Healpix_Map<double>& map);

    /// Calculate the power spectrum for a given map.
    /// \param mapName The name of the fits file containing the map.
    void calculate(const char* mapName);
    
    /// Retrieve the calculated power spectrum. Should be called after calculate.
    /// \return A constant reference to a map containing the power spectrum. The independent variable is l or the mid-points of bins if binning has been used.
    const std::map<double, double>& powerSpectrum() const { return ps_; }

    /// Retrieve the calculated pseudo power spectrum. Should be called after calculate.
    /// \return A constant reference to a vector containing the pseudo power spectrum. The index is l.
    const std::vector<double>& getPseudoSpectrum() const { return c_; }

    /// Retrieve the calculated pseude power spectrum for the mask.
    /// \return A constant reference to a vector containing the mask pseudo power spectrum. The index is l.
    const std::vector<double>& getMaskPseudoSpectrum() const { return w_; }

    /// Calculate and store the coupling kernel for a given mask.
    /// \param w A vector containing the weights of each pixel.
    /// \param lMax The maximum value of l.
    /// \param fileName The name of the file where the result should be stored.
    static void calculateCouplingKernel(const std::vector<double>& w, int lMax, const char* fileName);

    /// Calculate and store the coupling kernel for a given mask.
    /// \param mask The mask to be used. 
    /// \param lMax The maximum value of l.
    /// \param fileName The name of the file where the result should be stored.
    static void calculateCouplingKernel(const Healpix_Map<double>& mask, int lMax, const char* fileName);

    /// Calculate and store the coupling kernel for a given mask.
    /// \param maskName The name of the file containing the mask.
    /// \param lMax The maximum value of l.
    /// \param fileName The name of the file where the result should be stored.
    static void calculateCouplingKernel(const char* maskName, int lMax, const char* fileName);
    
private:
    void construct(const std::vector<int>* bins, int lMax);
    void calculate();
    void calculateCoupling();
    void calculateK();
    void calculatePS();
    
    double p(int b, int l) const;
    double q(int l, int b) const;
    
private:
    int lMax_;
    std::vector<int> bins_;
    std::vector<double> beam_;
    const Healpix_Map<double>* map_;
    Healpix_Map<double> mask_;
    Healpix_Map<double> maskedMap_;
    std::vector<double> c_, w_;
    std::string couplingKernelFileName_;
    
    std::vector<std::vector<double> > coupling_;
    LaGenMatDouble k_, kInv_;
    
    std::map<double, double> ps_;
};

#endif
