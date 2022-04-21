#include <string>
#include <vector>
#include <fstream>

std::size_t index(
  std::size_t x, 
  std::size_t y, 
  std::size_t z, 
  std::size_t depth, 
  std::size_t plane) {
  // 1D index corresponding to a flattened 3D variable. 
  return z + (y*depth) + (x*plane);
}

// void test_upper_dt(Options &options) {
//   float ka = options.k*options.a;
//   float k_1_a = options.k*(1 - options.a);
//   float max = (ka > k_1_a) ? ka : k_1_a;
//   float lambda = options.delta/(options.dx*options.dx);
//   float r_plus = options.k*(options.b+1)*(options.b+1)/4.0;
//   float upper_bound_dt = 1.0/(4*lambda + max + r_plus);
//   if (options.dt > upper_bound_dt) 
//     throw std::runtime_error(
//       "Forward Euler method is not stable, because dt ("+std::to_string(options.dt)+
//       ") > upper bound ("+std::to_string(upper_bound_dt)+")."
//     );
// }

int main (int argc, char** argv) {
  // Constants in the PDE computation
  const float my1 = 0.07;
  const float my2 = 0.3;
  const float delta = 5.0e-5;
  const float epsilon = 0.01;
  const float a = 0.1;
  const float b = 0.1;
  const float k = 8.0;
  const float dx = 0.000143;
  const float dt = 0.0001;

  // Mesh dimensions and number of iterations
  const std::size_t num_iterations = 10000;
  const std::size_t h = 50;
  const std::size_t w = 50;
  const std::size_t d = 50;

  // Padded mesh dimensions and respective volumes
  const std::size_t hp = h + 2;
  const std::size_t wp = w + 2;
  const std::size_t dp = d + 2;
  const std::size_t volume = h*w*d;
  const std::size_t volumep = hp*wp*dp;
  const std::size_t plane = w*d;
  const std::size_t planep = wp*dp;

  // Mesh variables
  std::vector<float> e_mesh(volumep);
  std::vector<float> e_temp(volume);
  std::vector<float> r_mesh(volume);
  float e_center;
  float r_center;

  // Other variables
  const float d_dx2 = delta/(dx*dx);
  std::size_t write_counter = 0; // Number of files that have been written
  std::size_t write_frequency = 150; // How frequent (number of iterations) to write
  std::size_t slice = d / 2; // Which slice index (along depth) to write

  // Initial e and r
  for (std::size_t x = 1; x <= h; ++x) {
    for (std::size_t y = 1; y <= w; ++y) {
      for (std::size_t z = 1; z <= d; ++z) {
        e_mesh[index(x,y,z,dp,planep)] = (y < w/2) ? 0.0 : 1.0;
        r_mesh[index(x-1,y-1,z-1,d,plane)] = (x < h/2) ? 1.0 : 0.0;
      }
    }
  }

  // Perform PDE model
  for (std::size_t t = 0; t < num_iterations; ++t) {
    /* 
    Write a slice from e_mesh and r_mesh.
    Files will be headerless csv-files of a slice in the mesh.
    Only write once per 'write_frequency' iteration.
    */
    if (t % write_frequency == 0) {
      std::string e_path = "./data/e" + std::to_string(write_counter) + ".csv";
      std::string r_path = "./data/r" + std::to_string(write_counter) + ".csv";
      write_counter++;
      std::ofstream e_file; 
      std::ofstream r_file; 
      e_file.open(e_path);
      r_file.open(r_path);

      // Save a slice
      for (std::size_t x = 1; x <= h; ++x) {
        for (std::size_t z = 1; z <= d; ++z) {
          // Print values
          e_file << e_mesh[index(x,slice,z,dp,planep)];
          r_file << r_mesh[index(x-1,slice-1,z-1,d,plane)];
          if (z < d) {
            e_file << ",";
            r_file << ",";
          }
        }
        e_file << "\n";
        r_file << "\n";
      }
      e_file.close();
      r_file.close();
    }

    /*
    Boundary condition: net-zero gradient.
    Enforced by copying inner neighbour surfaces to padded boundaries, 
    which mathematically is equivalent of zero flow out and in of the mesh.
    Since each dimension goes from 0, ..., h-1 (or w or d),
    the padded dimensions go from 0, ..., h+1 (or w or d).
    Immediate inner surfaces are at index 2 and h-1 (or w or d)
    The destiation index for those surfaces are at index 0 and h+1 (or w or d)
    */
    for (std::size_t x = 1; x <= h; ++x) {
      for (std::size_t y = 1; y <= w; ++y) {
        e_mesh[index(x,y,0,dp,planep)] = e_mesh[index(x,y,2,dp,planep)]; // front surface
        e_mesh[index(x,y,d+1,dp,planep)] = e_mesh[index(x,y,d-1,dp,planep)]; // back surface
      }
    }
    for (std::size_t x = 1; x <= h; ++x) {
      for (std::size_t z = 1; z <= d; ++z) {
        e_mesh[index(x,0,z,dp,planep)] = e_mesh[index(x,2,z,dp,planep)]; // left surface
        e_mesh[index(x,w+1,z,dp,planep)] = e_mesh[index(x,w-1,z,dp,planep)]; // right surface
      }
    }
    for (std::size_t y = 1; y <= w; ++y) {
      for (std::size_t z = 1; z <= d; ++z) {
        e_mesh[index(0,y,z,dp,planep)] = e_mesh[index(2,y,z,dp,planep)]; // top surface
        e_mesh[index(h+1,y,z,dp,planep)] = e_mesh[index(h-1,y,z,dp,planep)]; // bottom surface
      }
    }

    /*
    PDE computation by sliding stencil over inner volume 
    Note: inner volume of padded mesh corresponds to full volume of unpadded mesh, 
    and when calling the index function, correct width and depth must be used:
    * Variables w and d are the correct width and depth for r_mesh and e_temp.
    * Variables wp and dp are the correct width and depth for e_mesh.
    */
    for (std::size_t x = 1; x <= h; ++x) {
      for (std::size_t y = 1; y <= w; ++y) {
        for (std::size_t z = 1; z <= d; ++z) {
          e_center = e_mesh[index(x,y,z,dp,planep)]; // reusable variable
          r_center = r_mesh[index(x-1,y-1,z-1,d,plane)]; // reusable variable

          // New e_center (stored in e_temp)
          e_temp[index(x-1,y-1,z-1,d,plane)] = e_center + dt*(
            d_dx2*(-6*e_center + 
              e_mesh[index(x-1,y,z,dp,planep)] + e_mesh[index(x+1,y,z,dp,planep)] +
              e_mesh[index(x,y-1,z,dp,planep)] + e_mesh[index(x,y+1,z,dp,planep)] +
              e_mesh[index(x,y,z-1,dp,planep)] + e_mesh[index(x,y,z+1,dp,planep)]
            ) 
            - k*e_center*(e_center - a)*(e_center - 1) - e_center*r_center
          );

          // New r_center
          r_mesh[index(x-1,y-1,z-1,d,plane)] = r_center + dt*(
            -(epsilon + my1*r_center/(my2 + e_center))*
            (r_center + k*e_center*(e_center - b - 1))
          );
        }
      }
    }

    // Re-update e_mesh from e_temp
    for (std::size_t x = 1; x <= h; ++x) {
      for (std::size_t y = 1; y <= w; ++y) {
        for (std::size_t z = 1; z <= d; ++z) {
          e_mesh[index(x,y,z,dp,planep)] = e_temp[index(x-1,y-1,z-1,d,plane)];
        }
      }
    }
  }

  return EXIT_SUCCESS;
}
