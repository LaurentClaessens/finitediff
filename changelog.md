# Finitediff's changelog

## July 6, 2017 : the pointer returned by an exception's `what`

### Problem

`cppcheck` said (with some truth) that the method `what()` of an exception should not do
```
std::string text="foo"
return text.c_str()
```
because the returned pointer points to `text` which is deleted when exiting the function.

### Solution

Now a typical exception is structured as the following.

```
class NegativeMatrixElementNumberException : public std::exception
{
    private :
        std::string _msg;

        std::string message(const int n)
        {
            std::string s_num=std::to_string(n);
            return "Trying to access line or column with negative number : "+s_num;
        };

    public: 
        explicit NegativeMatrixElementNumberException(const int n): 
            _msg(message(n))
    {}
        virtual const char* what() const throw()
        {
            return _msg.c_str();
        }
};
```

- The lifetime of the object `_msg` is the lifetime of the exception object. So the pointer will be valid whenever one need him.
- When the exception object is destroyed, the string `_msg` is destroyed and there is no leak.
- The string is build (via the private function `message`) in the initiation list of the constructor because `what` has to be `const`.
