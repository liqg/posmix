#include "sam.h"
#include "faidx.h"
#include "interval_tree.hpp"  // 假设我们已经定义了 IntervalTree 类
#include <iostream>
#include <unordered_map>
#include <map>
#include <sstream>
#include <fstream>  // 包含 ifstream 的头文件
#include <iostream> // 包含 cout 和 cerr 等输入输出流的头文件
#include <string>
#include <vector>
#include <tuple>
#include <set>
#include <random>
#include <algorithm>
#include <getopt.h>
#include <regex>
#include "dataframe.hpp"

template <typename T>
std::vector<T> random_sample(const std::vector<T>& v, size_t k, unsigned int seed = 0) {
    std::vector<T> result(v);
    std::mt19937_64 rng(seed);
    std::shuffle(result.begin(), result.end(), rng);
    if (k < result.size()) result.resize(k);
    return result;
}
#include <string>

std::string replaceTabs(const std::string& str) {
    std::string result = str;
    std::string tab = "\t";
    std::size_t pos = result.find("\\t");
    while (pos != std::string::npos) {
        result.replace(pos, 2, tab);
        pos = result.find("\\t", pos + tab.length());
    }
    return result;
}

struct Config {
    std::string version = "v1.2";
    std::string bamfname1 = "";
    std::string bamfname2 = "";
    std::string outbamfname = "";
    std::string posfile = "";
    int seed = 1;
    int thread = 4;
    DataFrame<std::string, int, double> posdata;
    DataFrame<int, int, double> posstat;
    std::map<std::string, IntervalTree> trees;
    std::string newrg = "@RG\tID:tumor\tSM:tumor";
    std::string newrgid = "tumor";
    bool keepdup = false;
}; 

void tidyhdr(sam_hdr_t *hdr, const Config &conf) {
    if (conf.newrg != "") {
        sam_hdr_remove_lines(hdr, "RG", NULL, NULL);
        sam_hdr_add_lines(hdr, conf.newrg.c_str(), conf.newrg.size());
    }
}

void tidyaln(bam1_t *b, const Config &conf) {
    if (conf.newrg != "") bam_aux_update_str(b, "RG", conf.newrgid.size(), conf.newrgid.c_str());
    if (!conf.keepdup) b->core.flag &= ~BAM_FDUP;
}

bool check_chromosomes(const sam_hdr_t* hdr1, const sam_hdr_t* hdr2) {
    if (hdr1->n_targets != hdr2->n_targets) {
        return false;
    }
    for (int i = 0; i < hdr1->n_targets; ++i) {
        if (std::string(hdr1->target_name[i]) != std::string(hdr2->target_name[i])) {
            return false;
        }
    }
    return true;
}

void ReadPosFile(Config& conf) {
    std::ifstream file(conf.posfile);
    if (file) {
        std::string line;
        int id = 0;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::string chrom;
            int pos;
            double freq;
            iss >> chrom >> pos >> freq;
            if (!iss.fail()) {
                conf.posdata.add_row(std::make_tuple(chrom, pos, freq));
                conf.trees[chrom].insert(pos, pos+1, id);
                id++;
            } else {
                std::cerr << "Error reading line: " << line << '\n';
                exit(1);
            }
        }
    } else {
        std::cerr << "Error opening file\n";
        exit(1);
    }
}

// sequential bam reader
class BamSeqReader {
public:
    const char *fname = NULL;
    samFile *fp = NULL;
    sam_hdr_t *hdr = NULL;
    
    BamSeqReader (const char *fn) {
        fname = fn;
        fp = sam_open(fname, "r");
        if (!fp) {
            fprintf(stderr, "Couldn't open \"%s\" : %s", fname, strerror(errno));
            exit(EXIT_FAILURE);
        }
        hdr = sam_hdr_read(fp);
        if (!hdr) {
            fprintf(stderr, "Couldn't read header from \"%s\" : %s",
                    fname, strerror(errno));
            exit(EXIT_FAILURE);
        }
    }
    ~BamSeqReader() {
        if (hdr) {
            sam_hdr_destroy(hdr);
            hdr = NULL;
        }
        if (fp) {
            sam_close(fp);
            fp = NULL;
        }
     }
    int readaln(bam1_t *b) {
        int ret = sam_read1(fp, hdr, b);
        return ret;
     }
};

class BamWriter {
public:
    const char *fname;
    samFile *fp = NULL;
    sam_hdr_t *hdr = NULL;

    BamWriter(const char *fn, sam_hdr_t *hdr_) {
        fname = fn;
        hdr = hdr_;
        fp = sam_open(fname, "wb");
        if (!fp) {
            fprintf(stderr, "Couldn't open \"%s\" : %s", fname, strerror(errno));
            exit(EXIT_FAILURE);
        }

        if (sam_hdr_write(fp, hdr) < 0) {
            fprintf(stderr, "Couldn't write hdr into \"%s\" : %s", fname, strerror(errno));
            exit(EXIT_FAILURE);
        }
    }

    int writealn(bam1_t *b) {
        int ret = sam_write1(fp, hdr, b); 
        if ( ret < 0) {
            fprintf(stderr, "Couldn't write alignment into \"%s\" : %s", fname, strerror(errno));
            exit(EXIT_FAILURE);
        }
        return ret;
    }
    ~BamWriter(){
        sam_close(fp);
    }
};

class BamRegReader {
public:    
    const char *fname = NULL;
    htsFile *fp = NULL;
    hts_idx_t *idx = NULL;
    sam_hdr_t *hdr = NULL;
    hts_itr_t *itr = NULL;
    
    BamRegReader(const char *fn) {
        fname = fn;
        fp = sam_open(fname, "r");
        if (!fp) {
            fprintf(stderr, "Couldn't open \"%s\" : %s", fname, strerror(errno));
            exit(EXIT_FAILURE);
        }
        hdr = sam_hdr_read(fp);
        if (!hdr) {
            fprintf(stderr, "Couldn't read header from \"%s\" : %s",
                    fname, strerror(errno));
            exit(EXIT_FAILURE);
        }
        idx = sam_index_load(fp, fname);
        if (idx == NULL) {
            fprintf(stderr, "[%s] fail to load index for %s\n", __func__, fname);
            exit(EXIT_FAILURE);
        }        
     }
     ~BamRegReader() {
        if (itr) {
            sam_itr_destroy(itr);
            itr = NULL;
        }
        if (hdr) {
            sam_hdr_destroy(hdr);
            hdr = NULL;
        }
        if (idx) {
            hts_idx_destroy(idx);
            idx = NULL;
        }
        if (fp) {
            sam_close(fp);
            fp = NULL;
        }
     }

     void locate(const char *reg) {
        if (itr) {
            sam_itr_destroy(itr);
        }
        if ((itr=sam_itr_querys(idx, hdr, reg)) == 0) {
            fprintf(stderr, "[E::%s] fail to parse region %s\n", __func__, reg);
            exit(EXIT_FAILURE);
        }
     }

     int readaln(bam1_t *b) {
        if (itr == NULL) return -1;
        int ret = sam_itr_next(fp, itr, b);
        if (ret < 0) {
            sam_itr_destroy(itr);
            itr=NULL;
        }
        return ret;
     }
};

static void print_usage()
{
    fprintf(stderr,
"\n"
"Usage: posmix -1 receptor.bam -2 donor.bam -p posfreq.txt -o outmix.bam \n"
"\n"
"Input options:\n"
"  -1, --bam1,  receptor bam file\n"
"  -2, --bam2,  donor bam file\n"
"  -p, --pos,   position frequency file. Each line contains three tab-separated columns: chromosome, coordinate, and mixed  frequency. The mixed frequency refers to the proportion of reads of bam1 that are replaced by those from bam2 at this position.\n"
"\n"
"Output options:\n"
"  -o, --out,   output bam file\n"
"  -d, --keepdup, preserve the duplicate flag, which is removed by default.\n"
"  -r, --newrg, By default, all RG tags are removed and replaced with a new tag line: \"@RG\\tID:tumor\\tSM:tumor\". If the new tag is an empty string, the original tags are preserved.\n"
"\n"
"Generic options:\n"
"  -s, --seed,  seed for random sampling, default: 1\n"
"  -t, --thread,  number of threads for each bam reading/writing, default: 4\n"
);
}

int main(int argc, char** argv) {
    const char *optstr = "1:2:p:o:s:r:t:dhv";
    struct option opts[] = {
        {"bam1", 1, NULL, '1'},
        {"bam2", 1, NULL, '2'},
        {"pos", 1, NULL, 'p'},
        {"out", 1, NULL, 'o'},
        {"seed", 1, NULL, 's'},
        {"keepdup", 0, NULL, 'd'},
        {"newrg", 1, NULL, 'r'},
        {"thread", 1, NULL, 't'},
        {"version", 0, NULL, 'v'},
        {"help", 0, NULL, 'h'},
        {0, 0, 0, 0},
    };
    int opt;
    Config conf;
    while((opt = getopt_long(argc, argv, optstr, opts, NULL)) != -1){
        switch(opt) {
            case 'v':
                std::cerr << conf.version << std::endl;
                exit(0);
            case '1':
                conf.bamfname1 = optarg;
                break;
            case '2':
                conf.bamfname2 = optarg;
                break;
            case 'p':
                conf.posfile = optarg;
                break;
            case 'o':
                conf.outbamfname = optarg;
                break;
            case 's':
                conf.seed = std::atoi(optarg);
                break;
            case 't':
                conf.thread = std::atoi(optarg);
                break;
            case 'd':
                conf.keepdup = true;
                break;
            case 'r':
                conf.newrg = replaceTabs(std::string(optarg));
                //for (int i=0; i<7; i++) std::cerr << conf.newrg[i] << std::endl; exit(1);
                break;
            case 'h':
                print_usage();
                exit(0);
                break;
            case '?':
                if(strchr(optstr, optopt) == NULL){
                    fprintf(stderr, "unknown option '-%c'\n", optopt);
                    exit(EXIT_FAILURE);
                } else {
                    fprintf(stderr, "option requires an argument '-%c'\n", optopt);
                    exit(EXIT_FAILURE);
                }
                return 1;
        }
    }    
    if (strncmp(conf.newrg.c_str(), "@RG\tID:", 7) != 0) {
        std::cerr << "Error: invalid RG line: " << conf.newrg << std::endl;
        print_usage();
        exit(EXIT_FAILURE);
    }
    if (conf.newrg != "") {
        std::size_t l = conf.newrg.find("\t", 7) - 7;
        conf.newrgid = conf.newrg.substr(7, l);
        if(conf.newrgid == "") {
            fprintf(stderr, "Error: empty ID is not allowed.\n");
            exit(EXIT_FAILURE);
        }
    }

    if (conf.bamfname1 == "") {
        fprintf(stderr, "Error: missing background bam file\n");
        print_usage();
        exit(EXIT_FAILURE);
    }
    if (conf.bamfname2 == "") {
        fprintf(stderr, "Error: missing spike-in bam file\n");
        print_usage();
        exit(EXIT_FAILURE);
    }
    if (conf.outbamfname == "") {
        fprintf(stderr, "Eroor: missing output bam file\n");
        print_usage();
        exit(EXIT_FAILURE);
    }
    if (conf.posfile == "") {
        fprintf(stderr, "Error: missing position file\n");
        print_usage();
        exit(EXIT_FAILURE);
    }
    ReadPosFile(conf);

    auto bam1 = BamRegReader(conf.bamfname1.c_str());
    hts_set_threads(bam1.fp, conf.thread); // 设置使用4个线程
    auto bam2 = BamRegReader(conf.bamfname2.c_str());
    hts_set_threads(bam2.fp, conf.thread); // 设置使用4个线程
    if (!check_chromosomes(bam1.hdr, bam2.hdr)) {
        fprintf(stderr, "Error: Chromosome names and ordering in the headers of two bam files must be identical.\n");
        exit(EXIT_FAILURE);
    }
    std::unordered_map<std::string, int8_t> allnames1; // key: qname, v=1: discard, v=0: keep
    std::unordered_map<std::string, int8_t> allnames2; // key: qname, v=1: keep, v=0: discard

    bam1_t* b = bam_init1();
    for (int i = 0; i < conf.posdata.nrow(); ++i) {
        auto row = conf.posdata.get_row(i);
        std::string chr = std::get<0>(row);
        int pos = std::get<1>(row);
        double freq = std::get<2>(row);
        std::string reg = chr + ":" + std::to_string(pos-1) + "-" + std::to_string(pos+1);

        bam1.locate(reg.c_str());
        auto posnames1 = std::unordered_map<std::string, char>();
        while(bam1.readaln(b) >= 0) {
            std::string k = bam_get_qname(b);
            if (allnames1.find(k) == allnames1.end()) {
                posnames1.insert(std::make_pair(k, 0)); // 没有读过
            } else {
                posnames1.insert(std::make_pair(k, 1)); // 读过
            }
        }

        bam2.locate(reg.c_str());
        auto posnames2 = std::unordered_map<std::string, char>();
        while(bam2.readaln(b) >= 0) {
            std::string k = bam_get_qname(b);
            if (allnames2.find(k) == allnames2.end()) {
                posnames2.insert(std::make_pair(k, 0));
            } else {
                posnames2.insert(std::make_pair(k, 1));
            }
        }
        
        std::vector<std::string> newnames1;
        int nsample1 = 0.99 + posnames1.size()*freq;
        int nsample2 = nsample1;
        for (const auto& pair : posnames1) {
            if (pair.second == 0) {
                newnames1.emplace_back(pair.first);
            } else {
                nsample1 -= allnames1[pair.first];
            }
        }
        std::vector<std::string> newnames2;
        for (const auto& pair : posnames2) {
            if (pair.second == 0) {
                newnames2.emplace_back(pair.first);
            } else {
                nsample2 -= allnames2[pair.first]; // nsample--, if already been sampled
            }
        }

        nsample1 = std::min(int(newnames1.size()), nsample1);
        nsample2 = std::min(int(newnames2.size()), nsample2);
        nsample1 = std::max(0, nsample1);
        nsample2 = std::max(0, nsample2);

        for (auto i: random_sample(newnames1, nsample1, conf.seed)) {
            posnames1[i] = 2;
        }
        for (const auto& pair : posnames1) {// v=0: first read and keep, v=1: already in allnames1, v=2: first read and discard
            if (pair.second == 0) allnames1.insert(std::make_pair(pair.first, 0));
            if (pair.second == 2) allnames1.insert(std::make_pair(pair.first, 1));
        }

        for (auto i: random_sample(newnames2, nsample2, conf.seed)) {
            posnames2[i] = 2; // sampled
        }
        for (const auto& pair : posnames2) {
            if (pair.second == 0) allnames2.insert(std::make_pair(pair.first, 0));
            if (pair.second == 2) allnames2.insert(std::make_pair(pair.first, 1));
        }
    }
    bam_destroy1(b);

    // shrink allnames1 allnames2
    // only keep v=1, reads to discard
    for (auto it = allnames1.begin(); it != allnames1.end();) {
        // std::cerr << it->first << std::endl;
        if (it->second == 0) {
            it = allnames1.erase(it);
        } else {
            ++it;
        }
    }
    // only keep v=1, reads to keep
    for (auto it = allnames2.begin(); it != allnames2.end();) {
        if (it->second == 0) {
            it = allnames2.erase(it);
        } else {
            ++it;
        }
    }
    // std::cerr << "all pos processed" << std::endl;

    // output
    sam_hdr_t *outhdr = sam_hdr_dup(bam1.hdr);
    tidyhdr(outhdr, conf);

    auto outbam = BamWriter(conf.outbamfname.c_str(), outhdr);
    hts_set_threads(outbam.fp, conf.thread); // 设置使用4个线程
    
    // 按基因组坐标排序输出bam
    auto bamseq1 = BamSeqReader(conf.bamfname1.c_str());
    hts_set_threads(bamseq1.fp, conf.thread); // 设置使用4个线程
    bam1_t *b1 = bam_init1();
    int notover1 = bamseq1.readaln(b1) >= 0 ? true : false;
    auto bamseq2 = BamSeqReader(conf.bamfname2.c_str());
    hts_set_threads(bamseq2.fp, conf.thread); // 设置使用4个线程
    bam1_t *b2 = bam_init1();
    int notover2 = bamseq2.readaln(b2) >= 0 ? true : false;
    #define OUTPUTB1 \
        do {\
            tidyaln(b1, conf); \
            std::string k(bam_get_qname(b1));\
            if (allnames1.find(k) == allnames1.end()) outbam.writealn(b1); \
            notover1 = bamseq1.readaln(b1) >= 0 ? true : false; \
        } while(0)
    #define OUTPUTB2 \
        do {\
            tidyaln(b2, conf); \
            std::string k(bam_get_qname(b2)); \
            if (allnames2.find(k) != allnames2.end()) outbam.writealn(b2); \
            notover2 = bamseq2.readaln(b2) >= 0 ? true : false; \
        } while(0)

    while(1) {
        if (!notover1 && ! notover2) break;
        if (notover1 && !notover2) { OUTPUTB1; continue;}
        if (!notover1 && notover2) { OUTPUTB2; continue;}
        // if ((b1->core.flag & BAM_FUNMAP) != 0) {std::cerr << b1->core.tid << "\t" << b1->core.pos << "\t" << bam_get_qname(b1) << std::endl;}
        if (b1->core.tid < b2->core.tid ) { OUTPUTB1; continue; }
        if (b1->core.tid > b2->core.tid ) { OUTPUTB2; continue; }
        if (b1->core.pos <= b2->core.pos) { OUTPUTB1; continue; } else { OUTPUTB2; continue; }
    }
    bam_destroy1(b1);
    bam_destroy1(b2);

    return 0;
}