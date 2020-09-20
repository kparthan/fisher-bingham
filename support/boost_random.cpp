#include <boost/random.hpp>
#include <iostream> 
#include <ctime> 

//typedef boost::uniform_int<> NumberDistribution; 
typedef boost::uniform_real<> UniformReal; 
typedef boost::mt19937 RandomNumberGenerator; 
typedef boost::variate_generator<RandomNumberGenerator&,UniformReal> Generator; 

class UniformRandomNumberGenerator
{
  private:
    Generator numberGenerator;

  public:
    UniformRandomNumberGenerator(Generator &ex) : numberGenerator(ex)
    {}

    /*UniformRandomNumberGenerator operator=(const UniformRandomNumberGenerator &source)
    {
      if (this != &source) {
        numberGenerator = source.numberGenerator;
      }
      return *this;
    }*/

    Generator get()
    {
      return numberGenerator;
    }

    long double operator()()
    {
      //Generator numberGenerator(generator, distribution); 
      return numberGenerator();
    }
};

int main(int c, char** argv) 
{ 
  const int rangeMin = 0; 
  const int rangeMax = 1; 
  /*typedef boost::uniform_real<> NumberDistribution; 
  typedef boost::mt19937 RandomNumberGenerator; 
  typedef boost::variate_generator<RandomNumberGenerator&, 
                                   NumberDistribution> Generator; 
 
  NumberDistribution distribution(rangeMin, rangeMax); 
  RandomNumberGenerator generator; 
  Generator numberGenerator(generator, distribution); 
  generator.seed(std::time(0)); // seed with the current time 
 
  std::cout << numberGenerator() << std::endl; */

  UniformReal uniform_distribution(0,1);
  RandomNumberGenerator generator;
  Generator numberGenerator(generator,uniform_distribution); 
  generator.seed(std::time(0)); // seed with the current time 
  UniformRandomNumberGenerator instance(numberGenerator);
  for (int i=0; i<10; i++)
    std::cout << instance() << std::endl; 
  Generator ex = instance.get();
  UniformRandomNumberGenerator *copy = &instance;
  for (int i=0; i<10; i++)
    std::cout << (*copy)() << std::endl; 
} // main 
