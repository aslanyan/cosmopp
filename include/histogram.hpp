#ifndef COSMO_CPP_HISTOGRAM_HPP
#define COSMO_CPP_HISTOGRAM_HPP

#include <vector>
#include <map>
#include <limits>
#include <cmath>

namespace Math
{

/// A histogram class.

/// This class can be used to create a histogram from given data. A histogram is a map which maps the starting point of the bin to the bin value.
template<typename T>
class Histogram
{
public:
    /// The variable type
    typedef T VariableType;

    /// The histogram type
    typedef std::map<VariableType, int> HistogramType;
    
public:
    /// Constructor.

    /// Constructs an empty histogram. To generate a histogram at least one data point must be provided.
    Histogram() : binned_(false) { min_ = std::numeric_limits<VariableType>::max(); max_ = std::numeric_limits<VariableType>::min(); }
    
    /// Destructor.
    
    /// Destructor, does nothing.
    ~Histogram() {}
    
    /// Adds a data element.
    /// \param element The data element to be added.
    void addData(VariableType element);

    /// Get the minimum of the data points so far.
    /// \return The minimum.
    VariableType min() const { return min_; }

    /// Get the maximum of the data points so far.
    /// \return The maximum.
    VariableType max() const { return max_; }
    
    /// Get the number of data points provided so far.
    /// \return The data size.
    int getDataSize() const { return data_.size(); }
    
    /// This function creates an actual histogram out of the data given. Calling without any arguments creates a histogram on the whole range of data with an optimal number of bins.
    /// \param min The minimum of the range to be used.
    /// \param max The maximum of the range to be used.
    /// \param nBins The number of bins to be used.
    void createHistogram(VariableType min = std::numeric_limits<VariableType>::min(), VariableType max = std::numeric_limits<VariableType>::max(), int nBins = 0);
    
    /// Returns the actual histogram. Note that createHistogram needs to be called before calling this function.
    /// \return The histogram. It maps the starting points of the bins to the bin values.
    const HistogramType& getHistogram() const { check(binned_, "Histogram must be created first."); return histogram_; }
    
private:
    std::vector<VariableType> data_;
    VariableType min_;
    VariableType max_;
    
    HistogramType histogram_;
    bool binned_;
};
    
template<typename T>
void Histogram<T>::addData(VariableType element)
{
    if(binned_)
    {
        binned_ = false;
        histogram_.clear();
    }
    
    data_.push_back(element);
    if(element < min_)
        min_ = element;
    
    if(element > max_)
        max_ = element;
}
    
template<typename T>
void Histogram<T>::createHistogram(VariableType min, VariableType max, int nBins)
{
    check(!binned_, "Already created.");
    check(histogram_.empty(), "");
    check(!data_.empty(), "No data input yet.");
    
    if(min == std::numeric_limits<VariableType>::min())
        min = min_;
    if(max == std::numeric_limits<VariableType>::max())
        max = max_;
    
    check(max > min, "");
    
    if(nBins == 0)
        nBins = (int)std::ceil(std::sqrt((double)data_.size()));
    
    check(nBins > 0, "The number of bins must be positive.");
    
    VariableType step = (max - min) / nBins;
    for(int i = 0; i < nBins; ++i)
        histogram_[min + i * step] = 0;
    
    for(int i = 0; i < data_.size(); ++i)
    {
        if(data_[i] < min || data_[i] > max)
            continue;
        
        typename HistogramType::iterator it = histogram_.upper_bound(data_[i]);
        check(it != histogram_.begin(), "");
        
        --it;
        ++(*it).second;
    }
    
    binned_ = true;
}
	
} //namespace Math

#endif
