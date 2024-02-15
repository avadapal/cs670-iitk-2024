#include <assert.h>
#include <bsd/stdlib.h>
#include "dpf.h"

constexpr uint64_t nitems = 1ULL << (10);

using namespace dpf;


int main(int argc, char * argv[])
{
	AES_KEY prgkey;

	uint64_t target_value = 2;
	uint64_t target_index = 12;
	auto [dpfkey0, dpfkey1] = dpf_key<uint64_t, __m128i, AES_KEY>::gen(prgkey, nitems, target_index, target_value);

	
	for(size_t j = 0; j < nitems; ++j)
	{
		uint64_t x0 = eval(dpfkey0, j);
		uint64_t x1 = eval(dpfkey1, j);
		
		uint64_t reconstruction = x0 - x1;
		std::cout << j <<  "--> (reconstruction) " << reconstruction << " " << -reconstruction << std::endl;
	}	

	return 0;
}

