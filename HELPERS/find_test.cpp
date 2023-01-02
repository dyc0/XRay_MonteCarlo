#include <algorithm>
#include <cstdlib>
#include <vector>
#include <ctime>
#include <iostream>
#include <chrono>

using namespace std;

int find_index(double val, vector<double>& vec)
{
    for (int i=0; i<vec.size(); i++)
        if (vec[i] > val) return i;
    return -1;
}

int main()
{
    std::srand(std::time(nullptr));

    cout << "Generating vector" << endl;
    vector<double> vrec;
    for (int i=0; i<1e7; i++)
        vrec.push_back(rand());

    cout << "Sorting vector" << endl;
    sort(vrec.begin(), vrec.end());
    cout << "Vector generated." << endl;

    double lim = 1.*RAND_MAX/2;
    cout << "Lim: " << lim << endl;
    auto is_greater = [lim](double x) { return x > lim; };

    auto start = chrono::high_resolution_clock::now();
    std::cout << find_index(lim, vrec) << endl;
    auto stop = chrono::high_resolution_clock::now();
    auto duration1 = chrono::duration_cast<chrono::microseconds>(stop - start);
    cout << "MY DURATION: " << duration1.count() << endl;

    start = chrono::high_resolution_clock::now();
    std::cout << find_if(vrec.begin(), vrec.end(), is_greater) - vrec.begin() << endl;
    stop = chrono::high_resolution_clock::now();
    auto duration2 = chrono::duration_cast<chrono::microseconds>(stop - start);
    cout << "STD DURATION: " << duration2.count() << endl;

    cout << "DURATION COEFFICIENT: " << 1.*duration1.count()/duration2.count() << endl;

    return 0;
}
