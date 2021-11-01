#include <iostream>
#include <algorithm>
#include <vector>
#include <map>


struct Pair{
    int a;
    int b;

    Pair(int aa, int bb)
    {
        a = aa;
        b = bb;
    }
};


int main()
{
    std::map<int, std::vector<Pair>> indeces;
    unsigned int N;
    int S;
    std::cin >> N;
    std::vector<int> nums(N, 0);
    for(auto i=0; i<N; ++i){
        std::cin >> nums[i];
    }
    std::cin >> S;
    int sum;
    sort(nums.begin(), nums.end());
    for(int i=0; i<N-1; ++i){
        for(int j=i+1; j<N; ++j){
            sum = nums[i] + nums[j];
            indeces[sum].push_back(Pair(i, j));
        }
    }
    std::vector<int> four(4, 0);
    std::vector<Pair> pairs2;
    for(const auto &[sum, pairs1]: indeces){
        if(indeces.count(S-sum) == 1){
            pairs2 = indeces[S-sum];
            for(int k=0; k<pairs2.size(); ++k){
                if(pairs2[k].a != pairs1[].a && pairs2[k].a != j
                    && pairs2[k].b != i && pairs2[k].b != j){
                    four[0] = nums[pairs2[k].a];
                    four[1] = nums[pairs2[k].b];
                    four[2] = nums[i];
                    four[3] = nums[j];
                    sort(four.begin(), four.end());
                    std::cout << four[0] << ' ' << four[1] << ' ' <<
                                 four[2] << ' ' << four[3] << '\n';
                }
            }
        }
    }









//    for(int i=0; i<N-1; ++i){
//        for(int j=i+1; j<N; ++j){
//            sum = nums[i] + nums[j];
//            if(indeces.count(S-sum) == 1){
//                p = indeces[S-sum];
//                for(int k=0; k<p.size(); ++k){
//                    if(p[k].a != i && p[k].a != j && p[k].b != i && p[k].b != j){
//                        four[0] = nums[p[k].a];
//                        four[1] = nums[p[k].b];
//                        four[2] = nums[i];
//                        four[3] = nums[j];
//                        sort(four.begin(), four.end());
//                        std::cout << four[0] << ' ' << four[1] << ' ' <<
//                                     four[2] << ' ' << four[3] << '\n';
//                    }
//                }
//            }
//        }
//    }

    return 0;
}
