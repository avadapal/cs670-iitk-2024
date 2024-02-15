 
#ifndef DPFPP_DPF_H__
#define DPFPP_DPF_H__

#include <type_traits>  // std::is_same<>
#include <limits>       // std::numeric_limits<>
#include <climits>      // CHAR_BIT
#include <cmath>        // std::log2, std::ceil, std::floor
#include <stdexcept>    // std::runtime_error
#include <array>        // std::array<>
#include <iostream>     // std::istream and std::ostream
#include <vector>       // std::vector<>
#include <memory>       // std::shared_ptr<>
#include <utility>      // std::move
#include <algorithm>    // std::copy

#include <bsd/stdlib.h> // arc4random_buf
#include <x86intrin.h>  // SSE and AVX intrinsics

#include "aes.h"
#include "block.h"
#include "prg.h"

#define L 0
#define R 1
 
namespace dpf
{

template<typename leaf_t = __m128i, typename node_t = __m128i, typename prgkey_t = AES_KEY>
struct dpf_key;

template<typename leaf_t, typename node_t, typename prgkey_t>
inline leaf_t eval(const dpf_key <leaf_t, node_t, prgkey_t> & dpfkey, const size_t input);

template<typename node_t, typename prgkey_t>
static inline void expand(const prgkey_t & prgkey, const block<node_t> & seed, block<node_t> s[2], uint8_t t[2])
{
	dpf::PRG(prgkey, clear_lsb(seed, 0b11), s, 2);
	t[L] = get_lsb(s[L]);
	s[L] = clear_lsb(s[L], 0b11);
	t[R] = get_lsb(s[R]);
	s[R] = clear_lsb(s[R], 0b11);
} // dpf::expand

template<typename node_t, typename prgkey_t>
static inline void traverse(const prgkey_t & prgkey, const block<node_t> & seed, const bool direction,
	const uint8_t cw_t, const block<node_t> & cw, const uint8_t prev_t,
	block<node_t> & s, uint8_t & t)
{
	dpf::PRG(prgkey, clear_lsb(seed, 0b11), &s, 1, direction);
	t = get_lsb(s) ^ (cw_t & prev_t);
	//s = clear_lsb(xor_if( , , !prev_t), 0b11);
	// TODO: Fill up the first two parameters in the above function.
} // dpf::traverse

template<typename node_t, typename prgkey_t, size_t nodes_per_leaf>
static inline void stretch_leaf(const prgkey_t & prgkey, const block<node_t> & seed, std::array<block<node_t>, nodes_per_leaf> & s)
{
	dpf::PRG(prgkey, clear_lsb(seed), &s, nodes_per_leaf);
} // dpf::stretch_leaf

template<typename leaf_t, typename node_t, typename prgkey_t>
struct dpf_key final
{
  public:
	static constexpr size_t bits_per_leaf = std::is_same<leaf_t, bool>::value ? 1 : sizeof(leaf_t) * CHAR_BIT;
	static constexpr bool is_packed = (sizeof(leaf_t) < sizeof(node_t));
	static constexpr size_t leaves_per_node = dpf_key::is_packed ? sizeof(node_t) * CHAR_BIT / bits_per_leaf : 1;
	static constexpr size_t nodes_per_leaf = dpf_key::is_packed ? 1 : std::ceil(static_cast<double>(bits_per_leaf) / (sizeof(node_t) * CHAR_BIT));
	static_assert(leaves_per_node * bits_per_leaf == sizeof(node_t) * CHAR_BIT
	    || nodes_per_leaf * sizeof(node_t) == sizeof(leaf_t));

	using block_t = block<node_t>;
	using finalizer_t = std::array<block_t, nodes_per_leaf>;
	typedef std::pair<finalizer_t, finalizer_t> (*finalizer_callback)(const prgkey_t &, const size_t, const leaf_t &, const block_t[2], const uint8_t[2]);

	inline static constexpr size_t depth(const size_t nitems) { return std::ceil(std::log2(std::ceil(static_cast<double>(nitems) / dpf_key::leaves_per_node))); }
	inline constexpr size_t depth() const { return dpf_key::depth(nitems); }

	inline static constexpr size_t input_bits(const size_t nitems) { return std::ceil(std::log2(nitems)); }
	inline constexpr size_t input_bits() const { return dpf_key::input_bits(nitems); }

	inline dpf_key(dpf_key &&) = default;
	inline dpf_key & operator=(dpf_key &&) = default;
	inline dpf_key(const dpf_key &) = default;
	inline dpf_key & operator=(const dpf_key &) = default;

	inline bool operator==(const dpf_key & rhs) const { return nitems == rhs.nitems && root == rhs.root && cw == rhs.cw && finalizer == rhs.finalizer; }
	inline bool operator!=(const dpf_key & rhs) const { return !(*this == rhs); }

	static auto default_make_finalizer(const prgkey_t & prgkey, const size_t target, const leaf_t & val, const block_t s[2], const uint8_t t[2])
	{

		std::cout << "default+make+finalizer" << std::endl;
		finalizer_t finalizer;

		finalizer_t stretched[2];
		stretch_leaf(prgkey, s[L], stretched[L]);
		stretch_leaf(prgkey, s[R], stretched[R]);

		if constexpr(dpf_key::is_packed)
		{
			std::cout << "it is indeed packed! " << std::endl;
			auto finalizer0 = reinterpret_cast<block<node_t> *>(&finalizer[0]);
			if constexpr(std::numeric_limits<leaf_t>::is_integer)
			{
				if constexpr(std::is_same<leaf_t, bool>::value)
				{
					*finalizer0 = val ? 1 : 0;
				}
				else
				{
					typedef typename std::make_unsigned_t<leaf_t> unsigned_leaf_t;
					*finalizer0 = static_cast<unsigned_leaf_t>(val);
				}
				finalizer0->shiftl(bits_per_leaf * (target % leaves_per_node));
			}
			else
			{
				*finalizer0 = val;
			}
		}
		else
		{
			std::cout << "---> " << val[0] << " " << val[1] << std::endl;
			std::memcpy(&finalizer, &val, sizeof(finalizer_t));
		}
		for (size_t j = 0; j < nodes_per_leaf; ++j)
		{
			finalizer[j] ^= stretched[L][j] ^ stretched[R][j];
		}
		return std::make_pair(finalizer, finalizer);
	} // dpf_key::default_make_finalizer

	static auto gen(const prgkey_t & prgkey, size_t nitems, size_t target, const leaf_t & val = 1, const finalizer_callback make_finalizer = default_make_finalizer)
	{
		if (nitems <= target)
		{
			throw std::runtime_error("target point out of range");
		}

		block_t root[2];
		arc4random_buf(root, sizeof(root));
		uint8_t t[2] = { get_lsb(root[0]), !t[0] };
		root[1] = set_lsb(root[1], t[1]);
		block_t s[2] = { root[0], root[1] };

		const size_t depth = dpf_key::depth(nitems);
		std::vector<block_t> cw;
		cw.reserve(depth);

		block_t s0[2], s1[2];
		uint8_t t0[2], t1[2];
		const size_t nbits = input_bits(nitems);
		for (size_t layer = 0; layer < depth; ++layer)
		{
			const uint8_t bit = (target >> (nbits - layer - 1)) & 1U;

			expand(prgkey, s[0], s0, t0);
			expand(prgkey, s[1], s1, t1);

			const uint8_t keep = (bit == 0) ? L : R, lose = 1 - keep;
			bool cwt[2] = {
			    cwt[L] = t0[L] ^ t1[L] ^ bit ^ 1,
			    cwt[R] = t0[R] ^ t1[R] ^ bit
			};
			auto nextcw = s0[lose] ^ s1[lose];

			s[L] = xor_if(s0[keep], nextcw, !t[L]);
			t[L] = t0[keep] ^ (t[L] & cwt[keep]);

			s[R] = xor_if(s1[keep], nextcw, !t[R]);
			t[R] = t1[keep] ^ (t[R] & cwt[keep]);

			cw.emplace_back(set_lsbs(nextcw, cwt));
		}
		auto [finalizer0, finalizer1] = make_finalizer(prgkey, target, val, s, t);

		return std::make_pair(
		    std::forward<dpf_key>(dpf_key(nitems, root[0], cw, finalizer0, prgkey)),
		    std::forward<dpf_key>(dpf_key(nitems, root[1], cw, finalizer1, prgkey)));
	} // dpf_key::gen

	inline leaf_t eval(const size_t input) const { return std::forward<leaf_t>(dpf::eval(*this, input)); }
	
	const size_t nitems;
	const block_t root;
	const std::vector<block_t> cw;
	const finalizer_t finalizer;
	const prgkey_t prgkey;

  private:
	dpf_key(size_t nitems_, const block_t & root_, const std::vector<block_t> cw_,
		const finalizer_t & finalizer_, const prgkey_t & prgkey_)
	  : nitems(nitems_),
	    root(root_),
	    cw(cw_),
	    finalizer(finalizer_),
	    prgkey(prgkey_) { }
}; // struct dpf::dpf_key

 

	template<typename leaf_t, typename node_t>
	inline leaf_t getword(const block<node_t> & S, const size_t input)
	{
		auto S_ = reinterpret_cast<const leaf_t *>(&S);
		if constexpr(sizeof(leaf_t) >= sizeof(node_t)) return *S_;

		return S_[input];
	} // dpf::getword

	template<>
	inline bool getword(const block<__m128i> & S, const size_t input)
	{
		const __m128i mask = bool128_mask[input / 64];
		__m128i vcmp = _mm_xor_si128(_mm_and_si128(S >> (input % 64), mask), mask);

		return static_cast<bool>(_mm_testz_si128(vcmp, vcmp));
	} // dpf::getword<__m128i,bool>

	template<>
	inline bool getword(const block<__m256i> & S, const size_t input)
	{
		const __m256i mask = bool256_mask[input / 64];
		__m256i vcmp = _mm256_xor_si256(_mm256_and_si256(S >> (input % 64), mask), mask);

		return static_cast<bool>(_mm256_testz_si256(vcmp, vcmp));
	} // dpf::getword<__m256i,bool>

	template<typename leaf_t, typename node_t, typename prgkey_t>
	inline void finalize(const prgkey_t & prgkey, std::array<block<node_t>, dpf_key<leaf_t, node_t, prgkey_t>::nodes_per_leaf> finalizer, leaf_t * output, block<node_t> * s, size_t nnodes, uint8_t * t)
	{
		auto output_ = reinterpret_cast<std::array<block<node_t>, dpf_key<leaf_t, node_t, prgkey_t>::nodes_per_leaf> *>(output);
		for (size_t i = 0; i < nnodes; ++i)
		{
			stretch_leaf(prgkey, s[i], output_[i]);
			for (size_t j = 0; j < dpf_key<leaf_t, node_t, prgkey_t>::nodes_per_leaf; ++j)
			{
				output_[i][j] = xor_if(output_[i][j], finalizer[j], t[i]);
			}
		}
	} // dpf::finalize

 

	template<typename leaf_t, typename node_t, typename prgkey_t>
	inline leaf_t eval(const dpf_key<leaf_t, node_t, prgkey_t> & dpfkey, const size_t input)
	{
		auto prgkey = dpfkey.prgkey;
		auto root = dpfkey.root;
		auto depth = dpfkey.depth();
		auto nbits = dpfkey.input_bits();

		block<node_t> S = root;
		uint8_t T = get_lsb(root, 0b01);

		for (size_t layer = 0; layer < depth; ++layer)
		{
			auto & cw = dpfkey.cw[layer];
			const uint8_t nextbit = (input >> (nbits-layer-1)) & 1;
			traverse(prgkey, S, nextbit, get_lsb(cw, nextbit ? 0b10 : 0b01), cw, T, S, T); 
		}

		std::array<node_t, dpf_key<leaf_t, node_t, prgkey_t>::nodes_per_leaf> final;
		finalize(prgkey, dpfkey.finalizer, &final, &S, 1, &T);

		if constexpr(dpfkey.is_packed)
		{
			auto S_ = reinterpret_cast<block<node_t> *>(&final);
			return std::forward<leaf_t>(getword<leaf_t>(*S_, input % dpfkey.leaves_per_node));
		}
		else
		{
			auto ret = reinterpret_cast<leaf_t *>(&final);
			return *ret;
		}
	} // dpf::eval

} // namespace dpf

#endif
