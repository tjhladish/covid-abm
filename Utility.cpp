#include "Utility.h"

namespace covid {
    namespace util {
        Fit* lin_reg(const std::vector<double> &x, const std::vector<double> &y) {
            assert( x.size() == y.size() );
            Fit* fit = new Fit();
            const int n = x.size();
            double sumx = 0.0;                        // sum of x
            double sumx2 = 0.0;                       // sum of x**2
            double sumxy = 0.0;                       // sum of x * y
            double sumy = 0.0;                        // sum of y
            double sumy2 = 0.0;                       // sum of y**2

            for (int i=0; i<n; i++)   {
                sumx  += x[i];
                sumx2 += pow(x[i],2);
                sumxy += x[i] * y[i];
                sumy  += y[i];
                sumy2 += pow(y[i],2);
            }

            double denom = n * sumx2 - pow(sumx,2);
            if (denom == 0) {
                // singular matrix. can't solve the problem.
                fit->m   = 0;
                fit->b   = 0;
                fit->rsq = 0;
                return fit;
            }

            fit->m = (n * sumxy  -  sumx * sumy) / denom;
            fit->b = (sumy * sumx2  -  sumx * sumxy) / denom;
            // compute correlation coeff
            fit->rsq = pow((sumxy - sumx * sumy / n) / sqrt((sumx2 - pow(sumx,2)/n) * (sumy2 - pow(sumy,2)/n)),2);

            return fit;
        }

        istream& safeGetline(std::istream& is, std::string& t) {
            t.clear();

            // The characters in the stream are read one-by-one using a std::streambuf.
            // That is faster than reading them one-by-one using the std::istream.
            // Code that uses streambuf this way must be guarded by a sentry object.
            // The sentry object performs various tasks,
            // such as thread synchronization and updating the stream state.

            std::istream::sentry se(is, true);
            std::streambuf* sb = is.rdbuf();

            for(;;) {
                int c = sb->sbumpc();
                switch (c) {
                  case '\n':
                      return is;
                  case '\r':
                      if(sb->sgetc() == '\n')
                          sb->sbumpc();
                      return is;
                  case std::streambuf::traits_type::eof():
                      // Also handle the case when the last line has no line ending
                      if(t.empty())
                          is.setstate(std::ios::eofbit);
                      return is;
                  default:
                      t += (char)c;
                }
            }
        }

//        template <typename Out>
//        void split(const std::string &s, char delim, Out result) {
//            std::istringstream iss(s);
//            std::string item;
//            while (getline(iss, item, delim)) {
//                *result++ = item;
//            }
//        }
//
//        std::vector<std::string> split(const std::string &s, char delim) {
//            std::vector<std::string> elems;
//            split(s, delim, std::back_inserter(elems));
//            return elems;
//        }

        // TODO -- the split function below, and the commented out ones above, do not correctly return an empty string
        // if the empty string element is the last one on the line
        vector<string> split(const string &s, char delim) {
            vector<string> tokens;
            istringstream iss(s);
            string token;
            while (getline(iss, token, delim)) {
                tokens.push_back(token);
            }
            return tokens;
        }
    }
}
