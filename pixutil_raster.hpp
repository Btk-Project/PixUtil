#if !defined(_BTK_PROJECT_PIXUTIL_RASTER_HPP_)
#define _BTK_PROJECT_PIXUTIL_RASTER_HPP_

#include <algorithm>
#include <cmath>
//Raster algorithm
namespace PixRaster{
    #ifndef PIXUTIL_POINT
    struct Point{
        int x,y;
    };
    #else
    using Point = PIXUTIL_POINT;
    #endif

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
    template<typename Fn>
    void DrawVerticalLine(int x,int y1,int y2,Fn &&shade){
        if(y1 > y2){
            std::swap(y1,y2);
        }
        for(int y = y1;y <= y2;y++){
            shade(x,y);
        }
    }
    template<typename Fn>
    void DrawHorizontalLine(int x1,int x2,int y,Fn &&shade){
        if(x1 > x2){
            std::swap(x1,x2);
        }
        for(int x = x1;x <= x2;x++){
            shade(x,y);
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

    //Algo From http://www.sunshine2k.de/coding/java/TriangleRasterization/TriangleRasterization.html
    /**
     * @brief Fill Bottom-Left triangle with three points
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
    void FillBottomFlatTriangle(int x1,int y1,int x2,int y2,int x3,int y3,Fn &&shade,LineFn &&line){
        //       |
        //     |  |
        //    |    |
        //   |      |
        // |         |
        //|===========|
        float invslope1 = float(x2 - x1) / float(y2 - y1);
        float invslope2 = float(x3 - x1) / float(y3 - y1);

        float curx1 = x1;
        float curx2 = x1;

        //Begin scanline
        for(int scanline_y = y1;scanline_y <= y2;scanline_y++){
            line(curx1,scanline_y,curx2,scanline_y,shade);
            curx1 += invslope1;
            curx2 += invslope2;
        }
    }
    /**
     * @brief Fill Top-Left triangle with three points
     * 
     */
    template<typename Fn,typename LineFn>
    void FillTopFlatTriangle(int x1,int y1,int x2,int y2,int x3,int y3,Fn &&shade,LineFn &&line){
        //|=========|
        //  |      |
        //   |    |
        //    |  |
        //     |
        float invslope1 = float(x3 - x1) / float(y3 - y1);
        float invslope2 = float(x3 - x2) / float(y3 - y2);

        float curx1 = x3;
        float curx2 = x3;

        //Begin scanline
        for(int scanline_y = y3;scanline_y > y1;scanline_y--){
            line(curx1,scanline_y,curx2,scanline_y,shade);
            curx1 -= invslope1;
            curx2 -= invslope2;
        }
    }
    /**
     * @brief Generic fill a triangle
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
    void FillTriangle(int x1,int y1,int x2,int y2,int x3,int y3,Fn &&shade,LineFn &&line){
        Point ps[3] = {{x1,y1},{x2,y2},{x3,y3}};
        //Sort
        //TODO: use our own sort
        std::sort(ps,ps + 3,[](const Point &a,const Point &b){
            return a.y < b.y;
        });
        //ps[0].y < ps[1].y < ps[2].y
        if(ps[1].y == ps[2].y){
            FillBottomFlatTriangle(ps[0].x,ps[0].y,ps[1].x,ps[1].y,ps[2].x,ps[2].y,shade,line);
        }
        else if(ps[0].y == ps[1].y){
            FillTopFlatTriangle(ps[0].x,ps[0].y,ps[1].x,ps[1].y,ps[2].x,ps[2].y,shade,line);
        }
        else{
            //Split to two triangles
            Point p4 = {
                ps[0].x + (ps[1].y - ps[0].y) * (ps[2].x - ps[0].x) / (ps[2].y - ps[0].y),
                ps[1].y
            };
            FillBottomFlatTriangle(ps[0].x,ps[0].y,ps[1].x,ps[1].y,p4.x,p4.y,shade,line);
            FillTopFlatTriangle(ps[1].x,ps[1].y,p4.x,p4.y,ps[2].x,ps[2].y,shade,line);
        }
    }
    /**
     * @brief Fill a triangle with default line function
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
    void FillTriangle(int x1,int y1,int x2,int y2,int x3,int y3,Fn &&shade){
        FillTriangle(x1,y1,x2,y2,x3,y3,shade,DrawLine<Fn>);
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
    // template<typename Fn>
    // auto AddRotate(float x,float y,Fn &&shade){
    //     double f_cos = std::cos(x);
    //     double f_sin = std::sin(y);
    //     return [=](int in_x,int in_y){
    //         //Call shade
    //         shade(
    //             in_x * std::cos(y) - in_y * std::sin(y),
    //             in_x * std::sin(y) + in_y * std::cos(y)
    //         );
    //     };
    // }
} // namespace PixRaster


#endif // _BTK_PROJECT_PIXUTIL_RASTER_HPP_
