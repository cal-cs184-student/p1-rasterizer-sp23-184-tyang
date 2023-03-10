#include "rasterizer.h"

using namespace std;

namespace CGL {

  RasterizerImp::RasterizerImp(PixelSampleMethod psm, LevelSampleMethod lsm,
    size_t width, size_t height,
    unsigned int sample_rate) {
    this->psm = psm;
    this->lsm = lsm;
    this->width = width;
    this->height = height;
    this->sample_rate = sample_rate;


    sample_buffer.resize(width * height * sample_rate, Color::White);
  }

  // Used by rasterize_point and rasterize_line
  void RasterizerImp::fill_pixel(size_t x, size_t y, Color c) {
    // TODO: Task 2: You might need to this function to fix points and lines (such as the black rectangle border in test4.svg)
    // NOTE: You are not required to implement proper supersampling for points and lines
    // It is sufficient to use the same color for all supersamples of a pixel for points and lines (not triangles)

    //original implementation
    //sample_buffer[y * width + x] = c;
    for (int k = 0; k < 3; ++k) {
      this->rgb_framebuffer_target[3 * (y * width + x) + k] = (&c.r)[k] * 255;
    }
  }

  // helper function to add a buffer layer for supersamples
  void RasterizerImp::fill_sample_buffer(size_t x, size_t y, size_t index, Color c) {
    sample_buffer[sample_rate * (y * width + x) + index] = c;
  }

  // Rasterize a point: simple example to help you start familiarizing
  // yourself with the starter code.
  //
  void RasterizerImp::rasterize_point(float x, float y, Color color) {
    // fill in the nearest pixel
    int sx = (int)floor(x);
    int sy = (int)floor(y);

    // check bounds
    if (sx < 0 || sx >= width) return;
    if (sy < 0 || sy >= height) return;

    for (int index = 0; index < sample_rate; index++) {
      fill_sample_buffer(sx, sy, index, color);
    }
    return;
  }

  // Rasterize a line.
  void RasterizerImp::rasterize_line(float x0, float y0,
    float x1, float y1,
    Color color) {
    if (x0 > x1) {
      swap(x0, x1); swap(y0, y1);
    }

    float pt[] = { x0,y0 };
    float m = (y1 - y0) / (x1 - x0);
    float dpt[] = { 1,m };
    int steep = abs(m) > 1;
    if (steep) {
      dpt[0] = x1 == x0 ? 0 : 1 / abs(m);
      dpt[1] = x1 == x0 ? (y1 - y0) / abs(y1 - y0) : m / abs(m);
    }

    while (floor(pt[0]) <= floor(x1) && abs(pt[1] - y0) <= abs(y1 - y0)) {
      rasterize_point(pt[0], pt[1], color);
      pt[0] += dpt[0]; pt[1] += dpt[1];
    }
  }

  // Rasterize a triangle.
  void RasterizerImp::rasterize_triangle(float x0, float y0,
    float x1, float y1,
    float x2, float y2,
    Color color) {
    // TODO: Task 1: Implement basic triangle rasterization here, no supersampling
       
    // Create three triangle vertices
    //TODO: shift the vertices of the triangle to align with the scaled supersample pixel
    //TODO: if we shift the points within the vector from (x,y) to (x + (sqrt(sampling_rate) - 1)/2, y + (sqrt(sampling_rate) - 1)/2)
    Vector3D p0(x0, y0, 0);
    Vector3D p1(x1, y1, 0);
    Vector3D p2(x2, y2, 0);

    // create lines
    Vector3D line0 = p1 - p0;
    Vector3D line1 = p2 - p1;
    Vector3D line2 = p0 - p2;

    // create norms 
    Vector3D n0(-line0[1], line0[0], 0);
    Vector3D n1(-line1[1], line1[0], 0);
    Vector3D n2(-line2[1], line2[0], 0);

    // compute bounding box size
    int minX = floor(min({x0, x1, x2}));
    int maxX = ceil(max({x0,x1,x2}));
    int minY = floor(min({y0,y1,y2}));
    int maxY = ceil(max({y0,y1,y2}));

    // implement above line function
    int sample_rate_root = sqrt(sample_rate);
    for (int x = minX; x < maxX; x++) {
      for (int y = minY; y < maxY; y++) {
        int index = 0;
        Vector3D p(x, y, 0);
        // compute x, y coordinates based on sampling rate
        for (float i = 0.5; i < sample_rate_root; i++) {
          p.x = float(x) + i / sample_rate_root;
          for (float j = 0.5; j < sample_rate_root; j++) {
            p.y = float(y) + j / sample_rate_root;
            if (((dot(p - p0, n0) >= 0) && (dot(p - p1, n1) >= 0) && (dot(p - p2, n2) >= 0)) 
            || ((dot(p - p0, n0) < 0) && (dot(p - p1, n1) < 0) && (dot(p - p2, n2) < 0))) {
              fill_sample_buffer(x, y, index, color);
            }
            index++;
          }
        }
      }
    }
  }

    // TODO: Task 2: Update to implement super-sampled rasterization
    // IDEA: first, using the sampling rate, create a supersample buffer with the dimension width * height * sample_rate
    //such that each supersample is stored in the buffer.  Then, once this is done, we will go through the process of
    //sampling each point and including the color in the supersample buffer
    //after going through that process, we will then downsample by taking the average of all the values in the range [x, x+sample_rate]
    //to which we can write into our sample buffer.
    //TODO: find out what RasterizerImp::set_framebuffer_target and void RasterizerImp::resolve_to_framebuffer() do


  void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0,
    float x1, float y1, Color c1,
    float x2, float y2, Color c2)
  {
    // TODO: Task 4: Rasterize the triangle, calculating barycentric coordinates and using them to interpolate vertex colors across the triangle
    // Hint: You can reuse code from rasterize_triangle
      Vector3D p0(x0, y0, 0);
      Vector3D p1(x1, y1, 0);
      Vector3D p2(x2, y2, 0);

      // create lines
      Vector3D line0 = p1 - p0;
      Vector3D line1 = p2 - p1;
      Vector3D line2 = p0 - p2;

      // create norms
      Vector3D n0(-line0[1], line0[0], 0);
      Vector3D n1(-line1[1], line1[0], 0);
      Vector3D n2(-line2[1], line2[0], 0);

      // compute bounding box sizes for optimization
      int minX = floor(min({x0, x1, x2}));
      int maxX = ceil(max({x0,x1,x2}));
      int minY = floor(min({y0,y1,y2}));
      int maxY = ceil(max({y0,y1,y2}));

      // implement above line function
      int sample_rate_root = sqrt(sample_rate);
      for (int x = minX; x < maxX; x++) {
          for (int y = minY; y < maxY; y++) {
              int index = 0;
              Vector3D p(x, y, 0);
              // compute x, y coordinates based on sampling rate
              for (float i = 0.5; i < sample_rate_root; i++) {
                  p.x = float(x) + i / sample_rate_root;
                  for (float j = 0.5; j < sample_rate_root; j++) {
                      p.y = float(y) + j / sample_rate_root;
                      if (((dot(p - p0, n0) >= 0) && (dot(p - p1, n1) >= 0) && (dot(p - p2, n2) >= 0))
                          || ((dot(p - p0, n0) < 0) && (dot(p - p1, n1) < 0) && (dot(p - p2, n2) < 0))) {
                          Color newColor = interpolateColor(x0, y0, c0, x1, y1, c1, x2, y2, c2, p.x, p.y);
                          fill_sample_buffer(x, y, index, newColor);
                      }
                      index++;
                  }
              }
          }
      }
  }

  Color RasterizerImp::interpolateColor(float xA, float yA, Color c0,
                                       float xB, float yB, Color c1,
                                       float xC, float yC, Color c2,
                                       float x, float y) {
      float alpha = (-(x - xB)*(yC-yB) + ((y -yB)*(xC-xB)))/(-(xA - xB)*(yC-yB) + ((yA -yB)*(xC-xB)));
      float beta =  (-(x - xC)*(yA-yC) + ((y -yC)*(xA-xC)))/(-(xB - xC)*(yA-yC) + ((yB -yC)*(xA-xC)));
      float gamma = 1 - alpha - beta;
      float r = alpha * c0.r + beta * c1.r + gamma * c2.r;
      float g = alpha * c0.g + beta * c1.g + gamma * c2.g;
      float b = alpha * c0.b + beta * c1.b + gamma * c2.b;
      return {r, g, b};
  }


  void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0, float v0,
    float x1, float y1, float u1, float v1,
    float x2, float y2, float u2, float v2,
    Texture& tex)
  {
      //create the SampleParams struct
      SampleParams sp;
      sp.lsm = lsm;    
      sp.psm = psm;

      Vector3D p0(x0, y0, 0);
      Vector3D p1(x1, y1, 0);
      Vector3D p2(x2, y2, 0);

      // create lines
      Vector3D line0 = p1 - p0;
      Vector3D line1 = p2 - p1;
      Vector3D line2 = p0 - p2;

      // create norms
      Vector3D n0(-line0[1], line0[0], 0);
      Vector3D n1(-line1[1], line1[0], 0);
      Vector3D n2(-line2[1], line2[0], 0);

      //generate the new vectors and point to be used for the matrix
      Vector3D xyVec1 = p0 - p1;
      Vector3D xyVec2 = p2 - p1;
      //creating the vectors in origin
      Vector3D uvVec1(u0 - u1, v0 - v1, 0);
      Vector3D uvVec2(u2 - u1, v2 - v1, 0);
      Vector3D uvOrigin(u1, v1, 1);
      //create the matrix:
      Matrix3x3 xyMatrix(xyVec1.x, xyVec1.y, 0, xyVec2.x, xyVec2.y, 0, p1.x, p1.y, 1);
      Matrix3x3 uvMatrix(uvVec1.x, uvVec1.y, 0, uvVec2.x, uvVec2.y, 0, uvOrigin.x, uvOrigin.y, 1);
      Matrix3x3 xyToUV = uvMatrix * xyMatrix.inv();

      //basic matrix implementation
      Matrix3x3 M(x0, x1, x2, y0, y1, y2, 1, 1, 1);
      M = M.inv();
      Vector3D u = Vector3D(u0, u1, u2);
      Vector3D v = Vector3D(v0, v1, v2);

      // implement above line function
      int minX = floor(min({x0, x1, x2}));
      int maxX = ceil(max({x0,x1,x2}));
      int minY = floor(min({y0,y1,y2}));
      int maxY = ceil(max({y0,y1,y2}));

      int sample_rate_root = sqrt(sample_rate);
      for (int x = minX; x < maxX; x++) {
          for (int y = minY; y < maxY; y++) {
              int index = 0;
              Vector3D p(x, y, 1);
              // compute x, y coordinates based on sampling rate
              for (float i = 0.5; i < sample_rate_root; i++) {
                  p.x = float(x) + i / sample_rate_root;
                  for (float j = 0.5; j < sample_rate_root; j++) {
                      p.y = float(y) + j / sample_rate_root;
                      if (((dot(p - p0, n0) >= 0) && (dot(p - p1, n1) >= 0) && (dot(p - p2, n2) >= 0)) 
                      || ((dot(p - p0, n0) < 0) && (dot(p - p1, n1) < 0) && (dot(p - p2, n2) < 0))) { 
                        Vector3D uvPt = M * p;
                        Vector3D uv_dx = M * (p + Vector3D(1, 0, 0));
                        Vector3D uv_dy = M * (p + Vector3D(0, 1, 0));

                        Vector2D p_uv = Vector2D(dot(uvPt, u), dot(uvPt, v));
                        Vector2D p_uv_dx = Vector2D(dot(uv_dx, u), dot(uv_dx, v));
                        Vector2D p_uv_dy = Vector2D(dot(uv_dy, u), dot(uv_dy, v));

                        sp.p_uv = p_uv;
                        sp.p_dx_uv = p_uv_dx - p_uv;
                        sp.p_dy_uv = p_uv_dy - p_uv;
                        
                        Color color = tex.sample(sp);
                        fill_sample_buffer(x, y, index, color);
                      }
                      index++;
                  }
              }
          }
      }
    // TODO: Task 5: Fill in the SampleParams struct and pass it to the tex.sample function.
    //layout:
//      1. Write function that computes the matrix that takes the points xn,yn and converts to u,v basis ???> this can be done by writing MX = U and solving for M = UX^-1
//          1. First set up the matrix (x,y,0,x,y,0,x,y,1) and the matrix (u,v,0,u,v,0,u,v,1) and multiply X^-1 on the RHS
//      2. Using the matrix, convert the vector to a vector in the u,v space
//      3. For all points within the triangle, use the calculated matrix to generate the vector corresponding to it???s location within the u,v space.
//      4. Sample the pixel at that location in the texel space and write into the sample buffer that pixel
    // TODO: Task 6: Set the correct barycentric differentials in the SampleParams struct.
    // Hint: You can reuse code from rasterize_triangle/rasterize_interpolated_color_triangle
  }

  void RasterizerImp::set_sample_rate(unsigned int rate) {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->sample_rate = rate;
    //we need to create a supersample buffer in which we will store the increased dimension image
    //increasing the sample buffer size such that it now accounts each subpixel as a pixel
    this->sample_buffer.resize(width * height * rate, Color::White);
  }


  void RasterizerImp::set_framebuffer_target(unsigned char* rgb_framebuffer,
    size_t width, size_t height)
  {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->width = width;
    this->height = height;
    this->rgb_framebuffer_target = rgb_framebuffer;

    this->sample_buffer.resize(width * height * sample_rate, Color::White);
  }


  void RasterizerImp::clear_buffers() {
    std::fill(rgb_framebuffer_target, rgb_framebuffer_target + 3 * width * height, 255);
    std::fill(sample_buffer.begin(), sample_buffer.end(), Color::White);
  }


  // This function is called at the end of rasterizing all elements of the
  // SVG file.  If you use a supersample buffer to rasterize SVG elements
  // for antialising, you could use this call to fill the target framebuffer
  // pixels from the supersample buffer data.
  //
  //TODO: need to rewrite logic here such that we treat each block of sqrt(sample_rate) sqrt(sample_rate) as a block and do the arithmetic to deal with it
  void RasterizerImp::resolve_to_framebuffer() {
    //TODO: Task 2: You will likely want to update this function for supersampling support
    //original implementation before task 2
    //for (int x = 0; x < width; ++x) {
    //      for (int y = 0; y < height; ++y) {
    //          Color col = sample_buffer[y * width + x];
    //          for (int k = 0; k < 3; ++k) {
    //            this->rgb_framebuffer_target[3 * (y * width + x) + k] = (&col.r)[k] * 255;
    //        }
    //      }
    //    }
    for (int x = 0; x < width; ++x) {
      for (int y = 0; y < height; ++y) {
        Color col = getAvg(x,y);
        fill_pixel(x, y, col);
      }
    }
    //cout << "Completed resolving to framebuffer\n";
  }

  //helper function to compute average color based on sampling
  Color RasterizerImp::getAvg(int start_x, int start_y) {
        float r = 0;
        float g = 0;
        float b = 0;
        Color temp_color = {r, g, b};
        for (int index = 0; index < sample_rate; index++) {
          temp_color += sample_buffer[sample_rate * (start_y * width + start_x) + index];
        }

        temp_color.r /= sample_rate;
        temp_color.g /= sample_rate;
        temp_color.b /= sample_rate;

        return temp_color;
    }
    //this function allows for a user, given the sample_buffer, the original image width, the start sampling position of the pixel, and the sampling rate
    //to grab all of the information contained within the supersample square (x, y) to (x + sqrt(sampling_rate), y + sqrt(sampling_rate)
    //and outputs it at as a color object
  Rasterizer::~Rasterizer() { }


}// CGL
