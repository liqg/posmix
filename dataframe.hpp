#include <iostream>
#include <vector>
#include <tuple>

template<typename... Args>
class DataFrame {
public:
    DataFrame() {}

    void add_row(std::tuple<Args...> row) {
        data_.push_back(row);
    }

    std::tuple<Args...> get_row(int index) const {
        return data_[index];
    }

    int nrow() const {
        return data_.size();
    }

private:
    std::vector<std::tuple<Args...>> data_;
};

// int main() {
//     // 创建 data.frame
//     std::tuple<int, std::string, double> row1{1, "foo", 0.1};
//     std::tuple<int, std::string, double> row2{2, "bar", 0.2};
//     std::tuple<int, std::string, double> row3{3, "baz", 0.3};
//     DataFrame<int, std::string, double> df;
//     df.add_row(row1);
//     df.add_row(row2);
//     df.add_row(row3);

//     // 输出 data.frame
//     std::cout << "nrow: " << df.nrow() << std::endl;
//     for (int i = 0; i < df.nrow(); ++i) {
//         auto row = df.get_row(i);
//         std::cout << std::get<0>(row) << " " << std::get<1>(row) << " " << std::get<2>(row) << std::endl;
//     }

//     return 0;
// }