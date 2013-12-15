#ifndef COSMO_PP_EXCEPTION_HANDLER_HPP
#define COSMO_PP_EXCEPTION_HANDLER_HPP

#include <exception>
#include <string>

/// Basic exception with text.

/// This exception class is inherited from the C++ standard exception class and can just contain a text message.
class StandardException : public std::exception
{
public:
    /// Destructor.

    /// Destructor, does nothing.
    ~StandardException() throw() {}

    /// Function to get the text message.
    /// \return The text message for the exception.
    virtual const char* what() const throw()
    {
        return s_.c_str();
    }
    
    /// Function to set the text message.
    /// \param s Sets the text message.
    void set(const std::string& s) {s_ = s;}
    
private:
    std::string s_;
};

#endif
