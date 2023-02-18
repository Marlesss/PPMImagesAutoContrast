#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <time.h>
#include <omp.h>
#include <algorithm>

using namespace std;

int threads;

class PnmFile {
    public:
    string tof, sof, mcof;
    int width, height, max_color;
    vector <int> bytes;

    PnmFile(string filePath) {
        ifstream in(filePath, ios_base::binary);
        if (in.is_open()) {
            stringstream ss;
            string line;
            getline(in, line);
            tof = line;
            getline(in, line);
            sof = line;
            ss.clear();
            ss.str(line);
            ss >> width >> height;
            getline(in, line);
            mcof = line;
            ss.clear();
            ss.str(line);
            ss >> max_color;
            bytes = {istreambuf_iterator<char>(in), istreambuf_iterator<char>()};
            for (int i = 0; i < bytes.size(); i++) {
                if (bytes[i] < 0) {
                    bytes[i] += 256;
                }
            }
            in.close();
            if (!(tof == "P5" || tof == "P6")) {
                throw "Can't handle such type of file";
            }
            if (width < 0 || height < 0) {
                throw "Incorrect width or height of image";
            }
            if (max_color != 255) {
                throw "Incorrect max color. Waited for 255";
            }
        } else {
            throw "File \"" + filePath + "\" not found";
        }
    }

    void printInfo() {
        cout << "Type of file is " << tof << '\n';
        cout << "Width is " << width << '\n';
        cout << "Height is " << height << '\n';
        cout << "Max color is " << max_color << '\n';
    }

    void write(string filePath) {
        ofstream outfile(filePath, ofstream::binary);
        {
            ofstream::sentry sentry(outfile);
            if (sentry)
            {
                for (char x : tof) {
                    outfile.put(x);
                }
                outfile.put('\n');
                for (char x : sof) {
                    outfile.put(x);
                }
                outfile.put('\n');
                for (char x : mcof) {
                    outfile.put(x);
                }
                outfile.put('\n');
                for (int t : bytes) {
                    char x = t;
                    outfile.put(x);
                }
            } else {
                throw "Can't write in file \"" + filePath + "\"";
            }
        }
        outfile.close();
    }

    void autoBrightness(double coef) {
        if (coef < 0 || coef >= 0.5) {
            throw "Incorrect coefficient";
        }
        int bright[256];

        for (int i = 0; i < 255; i++) {
            bright[i] = 0;
        }
        #pragma omp parallel num_threads(threads)
        {
            #pragma omp for schedule(static)
            for (int i = 0; i < width * height; i++) {
                bright[brightnessOfPixel(i)]++;
            }
        }
        int ignore = (double) width * height * coef;
        int cnt_left = 0, left, cnt_right = 0, right;
        for (int i = 0; i < 256; i++) {
            cnt_left += bright[i];
            if (cnt_left > ignore) {
                left = i;
                break;
            }
        }
        for (int i = 255; i >= 0; i--) {
            cnt_right += bright[i];
            if (cnt_right > ignore) {
                right = i;
                break;
            }
        }
        left = min(max(left, 0), 255);
        right = min(max(right, 0), 255);
        if (left > right) {
            return;
        }
        #pragma omp parallel num_threads(threads)
        {
            #pragma omp for schedule(static)
            for (int i = 0; i < width * height; i++) {
                int old_brop = brightnessOfPixel(i);
                int new_brop = stretch(old_brop, left, right, 0, 255);
                changeBrightness(i, old_brop, new_brop);
            }
        }
    }

    void normalize(double coef) {
        int ignore = (double) width * height * coef;
        vector <int> brofpix(width*height);

        #pragma omp parallel for schedule(static)
        for (int i = 0; i < width * height; i++) {
            brofpix[i] = brightnessOfPixel(i);
        }
        sort(brofpix.begin(), brofpix.end());

        int x = brofpix[ignore];
        #pragma omp parallel num_threads(threads)
        {
            #pragma omp for schedule(static)
            for (int i = 0; i < bytes.size(); i++) {
                bytes[i] = max(bytes[i] - x, 0);
            }
        }

        x = brofpix[brofpix.size() - ignore - 1];
        double f = (double) 255 / x;
        #pragma omp parallel num_threads(threads)
        {
            #pragma omp for schedule(static)
            for (int i = 0; i < bytes.size(); i++) {
                bytes[i] = min(max((double) bytes[i] * f, 0.0), 255.0);
            }
        }
    }

    private:
    int brightnessOfPixel(int i) {
        if (tof == "P5") {
            return bytes[i];
        } else if (tof == "P6") {
            return (double) 0.3 * bytes[i * 3] + 0.59 * bytes[i * 3 + 1] + 0.11 * bytes[i * 3 + 2];
        }
        throw "Unknown type of file";
    }

    int stretch(double x, double from1, double from2, double to1, double to2) {
            double accurate = (x - from1) / (from2 - from1) * (to2 - to1) + to1;
            int result = accurate;
            return result;
    }

    void changeBrightness(int i, int old_brop, int new_brop) {
        if (tof == "P5") {
            bytes[i] = new_brop;
            return;
        } else if (tof == "P6") {
            if (old_brop != 0) {
                double k = min(min(min((double) new_brop / old_brop,
                                       (double) 255 / bytes[i * 3]),
                                   (double) 255 / bytes[i * 3 + 1]),
                               (double) 255 / bytes[i * 3 + 2]);
                bytes[i * 3] = (double) bytes[i * 3] * k;
                bytes[i * 3 + 1] = (double) bytes[i * 3 + 1] * k;
                bytes[i * 3 + 2] = (double) bytes[i * 3 + 2] * k;
            }
            return;
        }
        throw "Unknown type of file";
    }
};

int main(int argc, char* argv[])
{
    string inputPath;
    string outputPath;
    double coef;
    if (argc != 5) {
        cerr << "Error: Incorrect count of arguments\n";
        exit(1);
    }
    stringstream ss;
    ss.str(argv[1]);
    ss >> threads;
    ss.clear();
    ss.str(argv[2]);
    ss >> inputPath;
    ss.clear();
    ss.str(argv[3]);
    ss >> outputPath;
    ss.clear();
    ss.str(argv[4]);
    ss >> coef;
    ss.clear();
    try {
        PnmFile file(inputPath);
        clock_t start = clock();
        file.normalize(coef);
//        file.autoBrightness(coef);
        clock_t end = clock();
        file.write(outputPath);
        double seconds = (double)(end - start) / CLOCKS_PER_SEC * 1000;
        cout << "Time (" << threads << " thread(s)): " << seconds << " ms\n";
    } catch (const char* exception) {
        cerr << "Error: " << exception << "\n";
    }
}
