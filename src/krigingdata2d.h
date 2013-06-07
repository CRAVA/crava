/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/


#ifndef KRIGINGDATA2D_H
#define KRIGINGDATA2D_H

#include <vector>

class KrigingData2D
{
public:
  KrigingData2D(int nData = 0);
  ~KrigingData2D(void);

  void                       addData(int   i,
                                     int   j,
                                     float value);
  void                       findMeanValues(void);

  const std::vector<float> & getData(void)         const { return data_                          ;}
  const std::vector<int>   & getIndexI(void)       const { return indexI_                        ;}
  const std::vector<int>   & getIndexJ(void)       const { return indexJ_                        ;}

  inline int                 getNumberOfData(void) const { return static_cast<int>(data_.size()) ;}
  inline float               getData(int k)        const { return data_[k]                       ;}
  inline int                 getIndexI(int k)      const { return indexI_[k]                     ;}
  inline int                 getIndexJ(int k)      const { return indexJ_[k]                     ;}

  void                       writeToFile(const std::string & name);

private:
  int                        gotBlock(int i, int j) const;

  std::vector<float>         data_;
  std::vector<int>           indexI_;
  std::vector<int>           indexJ_;
  std::vector<int>           count_;

};

#endif
