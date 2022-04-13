# PixUtil  

PixUtil is a collection of utilities for processing images.  
written in C++17 and using templates.

| librarys | file                                    |  description |
| ---      | ---                                     | ---          |
| PixUtil  | [pixutil.hpp](pixutil.hpp)              | provide a view for pixel|
| PixImage | [pixutil_image.hpp](pixutil_image.hpp)  | provide a load / save   |
| PixFilter| [pixutil_fliter.hpp](pixutil_filter.hpp)| provide algorithms for filter|
| PixRaster| [pixutil_raster.hpp](pixutil_raster.hpp)| provide algorithms for raster|

## TODO List

- [ ] Add boundy check in Rotate

## Example

```cpp
#include "pixutil_filter.hpp"
#include "pixutil.hpp"

using namespace PixUtil;

int main(){
    //First create a view on pixel data etc.
    RGBAView view(/*Input your data here*/);
    //Then do your algorithm on the view
    view[y][x] = Pixels...;
    view[y][x] = Color{...};
    
    Filter2D(dst_view,view,YourFilter{...});
}

```
