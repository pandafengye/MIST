/**
 * @file    parseCmdArgs.hpp
 * @brief   Functionality related to command line parsing for indexing and mapping
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef PARSE_CMD_HPP 
#define PARSE_CMD_HPP

#include <iostream>
#include <string>
#include <fstream>
#include <cassert>

//Own includes
#include "map/include/map_parameters.hpp"
#include "map/include/map_stats.hpp"
#include "map/include/commonFunc.hpp"

//External includes
#include "common/argvparser.hpp"

namespace skch
{

  /**
   * @brief           Initialize the command line argument parser 
   * @param[out] cmd  command line parser object
   */
  void initCmdParser(CommandLineProcessing::ArgvParser &cmd)
  {
    cmd.setIntroductoryDescription("-----------------\n\
fastANI is a fast alignment-free implementation for computing whole-genome Average Nucleotide Identity (ANI) between genomes\n\
-----------------\n\
Example usage: \n\
$ fastANI -q genome1.fa -r genome2.fa -o output.txt\n\
$ fastANI -q genome1.fa --rl genome_list.txt -o output.txt");

    cmd.setHelpOption("h", "help", "Print this help page");

    cmd.defineOption("ref", "reference genome (fasta/fastq)[.gz]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("ref","r");

    cmd.defineOption("refList", "a file containing list of reference genome files, one genome per line", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("refList","rl");

    cmd.defineOption("query", "query genome (fasta/fastq)[.gz]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("query","q");

    cmd.defineOption("queryList", "a file containing list of query genome files, one genome per line", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("queryList","ql");

    cmd.defineOption("kmer", "kmer size <= 16 [default : 16]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("kmer","k");

    cmd.defineOption("threads", "thread count for parallel execution [default : 1]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("threads","t");

    cmd.defineOption("fragLen", "fragment length [default : 3,000]", ArgvParser::OptionRequiresValue);

    cmd.defineOption("minFraction", "minimum fraction of genome that must be shared for trusting ANI. If reference and query genome size differ, smaller one among the two is considered. [default : 0.2]", ArgvParser::OptionRequiresValue);

    cmd.defineOption("visualize", "output mappings for visualization, can be enabled for single genome to single genome comparison only [disabled by default]");

    cmd.defineOption("matrix", "also output ANI values as lower triangular matrix (format inspired from phylip). If enabled, you should expect an output file with .matrix extension [disabled by default]");

    cmd.defineOption("output", "output file name", ArgvParser::OptionRequired | ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("output","o");

    cmd.defineOption("version", "Show version", ArgvParser::NoOptionAttribute);
    cmd.defineOptionAlternative("version", "v");
  }

  /**
   * @brief                   Parse the file which has list of reference or query files
   * @param[in]   fileToRead  File containing list of ref/query files 
   * @param[out]  fileList    List of files will be saved in this vector   
   */
  template <typename VEC>
    void parseFileList(std::string &fileToRead, VEC &fileList)
    {
      std::string line;

      std::ifstream in(fileToRead);

      if (in.fail())
      {
        std::cerr << "ERROR, skch::parseFileList, Could not open " << fileToRead << "\n";
        exit(1);
      }

      while (std::getline(in, line))
      {
        //trim whitespaces
        skch::CommonFunc::trim (line);

        if (line.length() > 0)        //avoid empty strings
          fileList.push_back(line);
      }
    }

  /**
   * @brief                     validate the reference and query file(s)
   * @param[in] querySequences  vector containing query file names
   * @param[in] refSequences    vector containing reference file names
   */
  template <typename VEC>
    void validateInputFiles(VEC &querySequences, VEC &refSequences)
    {
      if (querySequences.size() == 0 || refSequences.size() == 0)
      {
        std::cerr << "ERROR, skch::validateInputFiles, Count of query and ref genomes should be non-zero" << std::endl;
        exit(1);
      }

      //Open file one by one
      for(auto &e : querySequences)
      {
        std::ifstream in(e);

        if (in.fail())
        {
          std::cerr << "ERROR, skch::validateInputFiles, Could not open " << e << std::endl;
          exit(1);
        }
      }

      for(auto &e : refSequences)
      {
        std::ifstream in(e);

        if (in.fail())
        {
          std::cerr << "ERROR, skch::validateInputFiles, Could not open " << e << std::endl;
          exit(1);
        }
      }
    }

  /**
   * @brief                   Print the parsed cmd line options
   * @param[in]  parameters   parameters parsed from command line
   */
  void printCmdOptions(skch::Parameters &parameters)
  {
    std::cerr << ">>>>>>>>>>>>>>>>>>" << std::endl;
    std::cerr << "Reference = " << parameters.refSequences << std::endl;
    std::cerr << "Query = " << parameters.querySequences << std::endl;
    std::cerr << "Kmer size = " << parameters.kmerSize << std::endl;
    std::cerr << "Fragment length = " << parameters.minReadLength << std::endl;
    std::cerr << "Threads = " << parameters.threads << std::endl;
    std::cerr << "ANI output file = " << parameters.outFileName << std::endl;
    std::cerr << ">>>>>>>>>>>>>>>>>>" << std::endl;
  }

  /**
   * @brief                   Parse the cmd line options
   * @param[in]   cmd
   * @param[out]  parameters  sketch parameters are saved here
   */
  void parseandSave(int argc, char** argv, 
      CommandLineProcessing::ArgvParser &cmd, 
      skch::Parameters &parameters)
  {
    int result = cmd.parse(argc, argv);

    //Make sure we get the right command line args
    if (cmd.foundOption("version"))
    {
      std::cerr << "version 1.3\n\n";
      exit(0);
    }
    if (result != ArgvParser::NoParserError)
    {
      std::cerr << cmd.parseErrorDescription(result) << "\n";

      if (result == ArgvParser::ParserHelpRequested)
        exit(0);
      else
        exit(1);
    }
    else if (!cmd.foundOption("ref") && !cmd.foundOption("refList"))
    { 
      std::cerr << "Provide reference file (s)\n";
      exit(1);
    }
    else if (!cmd.foundOption("query") && !cmd.foundOption("queryList"))
    { 
      std::cerr << "Provide query file (s)\n";
      exit(1);
    }

    std::stringstream str;

    //Parse reference files
    if(cmd.foundOption("ref"))
    {
      std::string ref;

      str << cmd.optionValue("ref");
      str >> ref;

      parameters.refSequences.push_back(ref);
    }
    else //list of files
    {
      std::string listFile;

      str << cmd.optionValue("refList");
      str >> listFile;

      parseFileList(listFile, parameters.refSequences);
    }

    //Size of reference
    //parameters.referenceSize = skch::CommonFunc::getReferenceSize(parameters.refSequences); 

    //fix reference length to a typical bacterial genome length
    parameters.referenceSize = 5000000;
    str.clear();

    //Parse query files
    if(cmd.foundOption("query"))
    {
      std::string query;

      str << cmd.optionValue("query");
      str >> query;

      parameters.querySequences.push_back(query);
    }
    else //list of files
    {
      std::string listFile;

      str << cmd.optionValue("queryList");
      str >> listFile;

      parseFileList(listFile, parameters.querySequences);
    }
    
    str.clear();

    parameters.alphabetSize = 4;

    parameters.reportAll = true;

    if(cmd.foundOption("visualize"))
    {
      if(parameters.refSequences.size() == 1 && parameters.querySequences.size() == 1)
      {
        parameters.visualize = true;
      }
      else
      {
         parameters.visualize = false;
         std::cerr << "WARNING, skch::parseandSave, visualization is disabled. It is not supported if more than one pair of genomes are asked to compare \n";
      }
    }
    else
      parameters.visualize = false;


    if(cmd.foundOption("matrix"))
      parameters.matrixOutput = true;
    else
      parameters.matrixOutput = false;

    //Parse algorithm parameters
    if(cmd.foundOption("kmer"))
    {
      str << cmd.optionValue("kmer");
      str >> parameters.kmerSize;
      str.clear();
    }
    else
    {
      if(parameters.alphabetSize == 4)
        parameters.kmerSize = 16;
      else
        parameters.kmerSize = 5;
    }

    if(cmd.foundOption("threads"))
    {
      str << cmd.optionValue("threads");
      str >> parameters.threads;
      str.clear();
    }
    else
      parameters.threads = 1;

    parameters.p_value = 1e-03;

    if(cmd.foundOption("fragLen"))
    {
      str << cmd.optionValue("fragLen");
      str >> parameters.minReadLength;
      str.clear();
    }
    else
      parameters.minReadLength = 3000;

    if(cmd.foundOption("minFraction"))
    {
      str << cmd.optionValue("minFraction");
      str >> parameters.minFraction;
      str.clear();
      assert(parameters.minFraction >= 0.0 && parameters.minFraction <= 1.0);
    }
    else
      parameters.minFraction = 0.2;

    parameters.percentageIdentity = 80;

    /*
     * Compute window size for sketching
     */

    //Compute optimal window size
    parameters.windowSize = skch::Stat::recommendedWindowSize(parameters.p_value,
        parameters.kmerSize, parameters.alphabetSize,
        parameters.percentageIdentity,
        parameters.minReadLength, parameters.referenceSize);

    str << cmd.optionValue("output");
    str >> parameters.outFileName;
    str.clear();

    printCmdOptions(parameters);

    //Check if files are valid
    validateInputFiles(parameters.querySequences, parameters.refSequences);
  }
}


#endif
