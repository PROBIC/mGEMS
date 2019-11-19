# cxxargs
Yet another C++11 header-only command line parsing library. Heavily
inspired by the [cxxopts
library](https://github.com/jarro2783/cxxopts), which you might want
to consider using.

## Accepted input
cxxargs will recognize command line arguments and options supplied in
the following formats
```
--long-name
--long-name=value
--long-name value

-a
-ab
-abc value
```
where the last option will take the argument 'value', while -a and -b
do not.

Anything after -- will be parsed as positional arguments (vector of
strings in the supplied order). Unrecognized arguments before -- will
be ignored.

# Usage
Simply drop the header file in your includes and add it to your project
```
#include "cxxargs.hpp"
```
## Example
Create a cxxargs::Arguments object
```
cxxargs::Arguments args("cxxargs command line parser", "Usage:
myprogram -d 0.5483 --toggle --list=1,2,3,4 -- positional1, positional2");
```
and define some possible arguments
```
  args.add_argument<double>('d', "double", "This is a double with the
  default value 0.222.", 0.222);
  args.add_long_argument<bool>("toggle", "This is a boolean toggle with
  the default value false.");
  args.add_short_argument<std::vector<int>>('l', "This is a list of
  integers with no default values.");
```
All arguments with no defaults are required and must be supplied from the command line.

Parse the command line input
```
myprogram --double=0.5483 --toggle -l 1,2,3,4 -- positional1, positional2
```
by calling
```
args.parse(argc, argv);
```
which will set the value of the double to 0.5483, switch the toggle to
true, create a vector of integers containing the values { 1, 2, 3,
4 }, and parse the positional arguments positional1 and positional2
into a vector of strings { "positional1", "positional2" }.

args.help() can be called
```
std::cout << args.help() << std::endl;
```
to print the help message
```
Usage: myprogram -d 0.5483 --toggle --list=1,2,3,4 -- positional1, positional2"
-d --double	This is a double with the default value 0.222.
-t --toggle	This is a boolean toggle with the default value false.
-l --list	This is a list of integers with no default values.
```

Parsed values can be accessed using
```
args.value<double>("double");
```
the value() function will throw an exception if 1) the argument was
defined without a default value and a value was not supplied to
parse(), or 2) the argument "double" has not defined as a
possible argument.

## Extensions
cxxargs can easily be extended to interpret the command-line input as something
unconvential. For example, to open the argument '-f infile.txt' as a
std::shared_ptr to a std::ifstream for reading, start by defining a std::istream&
operator>> for this particular class in the cxxargs namespace
```
namespace cxxargs {
  std::istream& operator>> (std::istream &in, std::shared_ptr<std::ifstream> &t) {
    std::string filename;
    in >> filename;
    t = std::shared_ptr<std::ifstream>(new std::ifstream(filename));
    return in;
  }
}
```
and use the same syntax as before to add the file
into the list of possible arguments
```
args.add_argument<std::shared_ptr<std::ifstream>>('f', "infile", "This is an infile stream.");
```

# Requirements and dependencies
This is a header-only library with no dependencies, only requiring
that your compiler supports C++11.

# License
The source code from this project is subject to the terms of the
Mozilla Public License, v. 2.0. A copy of the MPL is supplied with the
project, or can be obtained at
[https://mozilla.org/MPL/2.0/](https://mozilla.org/MPL/2.0/).
