#if !defined(_BTK_PROJECT_PIXUTIL_RASTER_HPP_)
#define _BTK_PROJECT_PIXUTIL_RASTER_HPP_

#include <algorithm>
#include <cmath>
//Raster algorithm
namespace PixRaster{
    /**
     * @brief Draw a line
     * 
     * @tparam Fn The shader function
     * @param x1 
     * @param y1 
     * @param x2 
     * @param y2 
     * @param shade 
     */
    template<typename Fn>
    void DrawLine(int x1,int y1,int x2,int y2,Fn &&shade){
        int dx = std::abs(x2 - x1);
        int dy = std::abs(y2 - y1);
        int sx = x1 < x2 ? 1 : -1;
        int sy = y1 < y2 ? 1 : -1;
        int err = dx - dy;
        while(true){
            shade(x1,y1);
            if(x1 == x2 && y1 == y2)
                break;
            int e2 = 2 * err;
            if(e2 > -dy){
                err -= dy;
                x1 += sx;
            }
            if(e2 < dx){
                err += dx;
                y1 += sy;
            }
        }
    }
    /**
     * @brief Draw a circle
     * 
     * @tparam Fn 
     * @param x 
     * @param y 
     * @param r 
     * @param shade 
     */
    template<typename Fn>
    void DrawCircle(int x,int y,int r,Fn &&shade){
        int x0 = r;
        int y0 = 0;
        int err = 0;
        while(x0 >= y0){
            shade(x + x0,y + y0);
            shade(x + y0,y + x0);
            shade(x - y0,y + x0);
            shade(x - x0,y + y0);
            shade(x - x0,y - y0);
            shade(x - y0,y - x0);
            shade(x + y0,y - x0);
            shade(x + x0,y - y0);
            y0++;
            err += 1 + 2 * y0;
            if(2 * (err - x0) + 1 > 0){
                x0--;
                err -= 2 * (x0 + 1);
            }
        }
    }
    /**
     * @brief Draw Triangle with three points
     * 
     * @tparam Fn 
     * @tparam LineFn 
     * @param x1 
     * @param y1 
     * @param x2 
     * @param y2 
     * @param x3 
     * @param y3 
     * @param shade 
     * @param line 
     */
    template<typename Fn,typename LineFn>
    void DrawTriangle(int x1,int y1,int x2,int y2,int x3,int y3,Fn &&shade,LineFn &&line){
        line(x1,y1,x2,y2,shade);
        line(x2,y2,x3,y3,shade);
        line(x3,y3,x1,y1,shade);
    }
    /**
     * @brief Draw triangle with default line function
     * 
     * @tparam Fn 
     * @param x1 
     * @param y1 
     * @param x2 
     * @param y2 
     * @param x3 
     * @param y3 
     * @param shade 
     */
    template<typename Fn>
    void DrawTriangle(int x1,int y1,int x2,int y2,int x3,int y3,Fn &&shade){
        DrawTriangle(x1,y1,x2,y2,x3,y3,shade,DrawLine<Fn>);
    }

} // namespace PixRaster

//Some useful function template
namespace PixRaster{
    template<typename Fn>
    auto AddScissor(int x,int y,int w,int h,Fn &&shade){
        return [=](int in_x,int in_y){
            //Test the point is in the scissor
            if(in_x < x || in_x >= x + w || in_y < y || in_y >= y + h)
                return;
            //Call shade
            shade(in_x,in_y);
        };
    }
    template<typename Fn>
    auto AddTranslate(int x,int y,Fn &&shade){
        return [=](int in_x,int in_y){
            //Call shade
            shade(in_x + x,in_y + y);
        };
    }
    template<typename Fn>
    auto AddScale(float x_f,float y_f,Fn &&shade){
        return [=](int in_x,int in_y){
            //Call shade
            shade(in_x * x_f,in_y * y_f);
        };
    }
    template<typename Fn>
    auto AddRotate(int x,int y,Fn &&shade){
        return [=](int in_x,int in_y){
            //Call shade
            shade(
                x + in_x * std::cos(y) - in_y * std::sin(y),
                y + in_x * std::sin(y) + in_y * std::cos(y)
            );
        };
    }
} // namespace PixRaster


#endif // _BTK_PROJECT_PIXUTIL_RASTER_HPP_
