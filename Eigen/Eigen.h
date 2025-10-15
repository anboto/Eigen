// SPDX-License-Identifier: Apache-2.0
// Copyright 2021 - 2025, the Anboto author and contributors
#ifndef _Eigen_Eigen_h
#define _Eigen_Eigen_h

#define EIGEN_MATRIX_PLUGIN 	<Eigen/ToStringPlugin.h>
#define EIGEN_DENSEBASE_PLUGIN 	<Eigen/ToStringPlugin.h>
#define EIGEN_TENSOR_PLUGIN		<Eigen/ToStringPlugin.h>

#define EIGEN_MPL2_ONLY

#ifndef _DEBUG
#define EIGEN_NO_DEBUG
#else
#define EIGEN_INITIALIZE_MATRICES_BY_NAN 
#endif

#define eigen_assert(x) ASSERT(x)

//  Re-export NEON types and intrinsics from U++ into the global namespace before including Eigen headers
#ifdef CPU_NEON	
namespace {
    using Upp::float64x1_t;
    using Upp::float64x2_t;
    using Upp::float32x4_t;
    using Upp::float32x2_t;
    using Upp::float32x4_t;
    using Upp::int8x8_t;
    using Upp::int8x16_t;
    using Upp::uint8x8_t;
    using Upp::uint8x16_t;
    using Upp::int16x4_t;
    using Upp::int16x8_t;
    using Upp::uint16x4_t;
    using Upp::uint16x8_t;
    using Upp::int32x2_t;
    using Upp::int32x4_t;
    using Upp::uint32x2_t;
    using Upp::uint32x4_t;
    using Upp::int64x2_t;
    using Upp::uint64x2_t;
    using Upp::vdup_n_f32;
    using Upp::vdupq_n_f32;
    using Upp::vreinterpret_s32_s8;
    using Upp::vdup_n_s8;
    using Upp::vget_low_s32;
    using Upp::vdupq_n_s8;
    using Upp::vdup_n_u8;
    using Upp::vreinterpret_u32_u8;
    using Upp::vdup_n_u8;
    using Upp::vdupq_n_u8;
    using Upp::vdup_n_s16;
    using Upp::vdupq_n_s16;
    using Upp::vdup_n_u16;
    using Upp::vdupq_n_u16;
    using Upp::vdup_n_s32;
    using Upp::vdupq_n_s32;
    using Upp::vdup_n_u32;
    using Upp::vdupq_n_u32;
    using Upp::vdupq_n_s64;
    using Upp::vdupq_n_u64;
    using Upp::vdup_n_u32;
    using Upp::vreinterpret_f32_u32;
    using Upp::vreinterpretq_f32_u32;
    using Upp::vadd_f32;
	using Upp::vaddq_f32;
	using Upp::vadd_s8;
	using Upp::vaddq_s8;
	using Upp::vadd_u8;
	using Upp::vaddq_u8;
	using Upp::vadd_s16;
	using Upp::vadd_u16;
	using Upp::vaddq_s16;
	using Upp::vaddq_u16;
	using Upp::vadd_s32;
	using Upp::vaddq_s32;
	using Upp::vadd_u32;
	using Upp::vaddq_u32;
	using Upp::vaddq_s64;
	using Upp::vaddq_u64;
	using Upp::vsub_f32;
	using Upp::vsubq_f32;
	using Upp::vsub_s8;
	using Upp::vsubq_s8;
	using Upp::vsub_u8;
	using Upp::vsubq_u8;
	using Upp::vsub_s16;
	using Upp::vsubq_s16;
	using Upp::vsub_u16;
	using Upp::vsubq_u16;
	using Upp::vsub_s32;
	using Upp::vsubq_s32;
	using Upp::vsub_u32;
	using Upp::vsubq_u32;
	using Upp::vsub_s64;
	using Upp::vsubq_s64;
	using Upp::vsub_u64;
	using Upp::vsubq_u64;
	using Upp::vneg_f32;
	using Upp::vnegq_f32;
	using Upp::vneg_s8;
	using Upp::vnegq_s8;
	using Upp::vneg_s16;
	using Upp::vnegq_s16;
	using Upp::vneg_s32;
	using Upp::vnegq_s32;
	using Upp::vnegq_s64;
	using Upp::vmul_f32;
	using Upp::vmulq_f32;
	using Upp::vmul_s8;
	using Upp::vmulq_s8;
	using Upp::vmul_u8;
	using Upp::vmulq_u8;
	using Upp::vmul_s16;
	using Upp::vmulq_s16;
	using Upp::vmul_u16;
	using Upp::vmulq_u16;
	using Upp::vmul_s32;
	using Upp::vmulq_s32;
	using Upp::vmul_u32;
	using Upp::vmulq_u32;
	using Upp::vdup_n_s64;
	using Upp::vcombine_s64;
	using Upp::vdup_n_u64;
	using Upp::vcombine_u64;
	using Upp::vreinterpret_s8_u32;
	using Upp::vreinterpret_u8_u32;
	using Upp::vreinterpret_s8_s32;
	using Upp::vreinterpret_u8_s32;
	using Upp::vfmaq_f32;
	using Upp::vfma_f32;
	using Upp::vfmsq_f32;
	using Upp::vfms_f32;
	using Upp::vmla_s8;
	using Upp::vmlaq_s8;
	using Upp::vmla_u8;
	using Upp::vmlaq_u8;
	using Upp::vmla_s16;
	using Upp::vmlaq_s16;
	using Upp::vmla_u16;
	using Upp::vmlaq_u16;
	using Upp::vmla_s32;
	using Upp::vmlaq_s32;
	using Upp::vmla_u32;
	using Upp::vmlaq_u32;
	using Upp::vabd_f32;
	using Upp::vabdq_f32;
	using Upp::vabd_s8;
	using Upp::vabdq_s8;
	using Upp::vabdq_s8;
	using Upp::vreinterpret_u32_u8;
	using Upp::vreinterpret_u8_u32;
	using Upp::vabd_u8;
	using Upp::vabdq_u8;
	using Upp::vabd_s16;
	using Upp::vabdq_s16;
	using Upp::vabd_u16;
	using Upp::vabdq_u16;
	using Upp::vabd_s32;
	using Upp::vabdq_s32;
	using Upp::vabd_u32;
	using Upp::vabdq_u32;
	using Upp::vmin_f32;
	using Upp::vminq_f32;
	using Upp::vminnmq_f32;
	using Upp::vminnm_f32;
	using Upp::vreinterpret_s32_s8;
	using Upp::vreinterpret_s8_s32;
	using Upp::vmin_s8;
	using Upp::vminq_s8;
	using Upp::vreinterpret_u32_u8;
	using Upp::vreinterpret_u8_u32;
	using Upp::vmin_u8;
	using Upp::vminq_u8;
	using Upp::vmin_s16;
	using Upp::vminq_s16;
	using Upp::vmin_u16;
	using Upp::vminq_u16;
	using Upp::vmin_s32;
	using Upp::vminq_s32;
	using Upp::vmin_u32;
	using Upp::vminq_u32;
	using Upp::vmax_f32;
	using Upp::vmaxq_f32;
	using Upp::vmaxnmq_f32;
	using Upp::vmaxnm_f32;
	using Upp::vmax_s8;
	using Upp::vmaxq_s8;
	using Upp::vmax_u8;
	using Upp::vmaxq_u8;
	using Upp::vmax_s16;
	using Upp::vmaxq_s16;
	using Upp::vmax_u16;
	using Upp::vmaxq_u16;
	using Upp::vmax_s32;
	using Upp::vmaxq_s32;
	using Upp::vmax_u32;
	using Upp::vmaxq_u32;
	using Upp::vcle_f32;
	using Upp::vcleq_f32;
	using Upp::vreinterpret_s32_u8;
	using Upp::vreinterpret_s8_u32;
	using Upp::vreinterpret_s8_u8;
	using Upp::vreinterpretq_s8_u8;
	using Upp::vcle_s8;
	using Upp::vcleq_s8;
	using Upp::vcle_u8;
	using Upp::vcleq_u8;
	using Upp::vcle_s16;
	using Upp::vcleq_s16;
	using Upp::vcle_u16;
	using Upp::vcleq_u16;
	using Upp::vcle_s32;
	using Upp::vcleq_s32;
	using Upp::vcle_u32;
	using Upp::vcleq_u32;
	using Upp::vreinterpret_s16_u16;
	using Upp::vreinterpretq_s16_u16;
	using Upp::vreinterpret_s32_u32;
	using Upp::vreinterpretq_s32_u32;
	using Upp::vcleq_s64;
	using Upp::vcleq_u64;
	using Upp::vreinterpretq_s64_u64;
	using Upp::vclt_f32;
	using Upp::vcltq_f32;
	using Upp::vclt_s8;
	using Upp::vcltq_s8;
	using Upp::vclt_u8;
	using Upp::vcltq_u8;
	using Upp::vclt_s16;
	using Upp::vcltq_s16;
	using Upp::vclt_u16;
	using Upp::vcltq_u16;
	using Upp::vclt_s32;
	using Upp::vcltq_s32;
	using Upp::vclt_u32;
	using Upp::vcltq_u32;
	using Upp::vclt_s64;
	using Upp::vcltq_s64;
	using Upp::vclt_u64;
	using Upp::vcltq_u64;
	using Upp::vceq_f32;
	using Upp::vceqq_f32;
	using Upp::vceq_s8;
	using Upp::vceqq_s8;
	using Upp::vceq_u8;
	using Upp::vceqq_u8;
	using Upp::vceq_s16;
	using Upp::vceqq_s16;
	using Upp::vceq_u16;
	using Upp::vceqq_u16;
	using Upp::vceq_s32;
	using Upp::vceqq_s32;
	using Upp::vceq_u32;
	using Upp::vceqq_u32;
	using Upp::vceqq_s64;
	using Upp::vceqq_u64;
	using Upp::vcge_f32;
	using Upp::vmvn_u32;
	using Upp::vcgeq_f32;
	using Upp::vmvnq_u32;
	using Upp::vreinterpret_u32_f32;
	using Upp::vreinterpretq_u32_f32;
	using Upp::vand_u32;
	using Upp::vandq_u32;
	using Upp::vand_s8;
	using Upp::vandq_s8;
	using Upp::vand_u8;
	using Upp::vandq_u8;
	using Upp::vand_s16;
	using Upp::vandq_s16;
	using Upp::vand_u16;
	using Upp::vandq_u16;
	using Upp::vand_s32;
	using Upp::vandq_s32;
	using Upp::vand_u32;
	using Upp::vandq_u32;
	using Upp::vandq_s64;
	using Upp::vandq_u64;
	using Upp::vorr_u32;
	using Upp::vreinterpret_u32_f32;
	using Upp::vreinterpret_f32_u32;
	using Upp::vorrq_u32;
	using Upp::vreinterpretq_u32_f32;
	using Upp::vreinterpretq_f32_u32;
	using Upp::vorr_s8;
	using Upp::vorrq_s8;
	using Upp::vorr_u8;
	using Upp::vorrq_u8;
	using Upp::vorr_s16;
	using Upp::vorrq_s16;
	using Upp::vorr_u16;
	using Upp::vorrq_u16;
	using Upp::vorr_s32;
	using Upp::vorrq_s32;
	using Upp::vorr_u32;
	using Upp::vorrq_u32;
	using Upp::vorrq_s64;
	using Upp::vorrq_u64;
	using Upp::veor_u32;
	using Upp::veorq_u32;
	using Upp::veor_s8;
	using Upp::veorq_s8;
	using Upp::veor_u8;
	using Upp::veorq_u8;
	using Upp::veor_s16;
	using Upp::veorq_s16;
	using Upp::veor_u16;
	using Upp::veorq_u16;
	using Upp::veor_s32;
	using Upp::veorq_s32;
	using Upp::veor_u32;
	using Upp::veorq_u32;
	using Upp::veorq_s64;
	using Upp::veorq_u64;
	using Upp::vbic_u32;
	using Upp::vbicq_u32;
	using Upp::vbic_s8;
	using Upp::vbicq_s8;
	using Upp::vbic_u8;
	using Upp::vbicq_u8;
	using Upp::vbic_s16;
	using Upp::vbicq_s16;
	using Upp::vbic_u16;
	using Upp::vbicq_u16;
	using Upp::vbic_s32;
	using Upp::vbicq_s32;
	using Upp::vbic_u32;
	using Upp::vbicq_u32;
	using Upp::vbicq_s64;
	using Upp::vbicq_u64;
	using Upp::vreinterpret_s32_s8;
	using Upp::vreinterpret_s8_s32;
	using Upp::vdup_n_s32;
	using Upp::vreinterpret_s32_u8;
	using Upp::vreinterpret_u8_s32;
	using Upp::vreinterpret_s8_u8;
	using Upp::vreinterpretq_s8_u8;
	using Upp::vreinterpret_u16_s16;
	using Upp::vreinterpretq_u16_s16;
	using Upp::vreinterpret_s32_u32;
	using Upp::vreinterpret_u32_s32;
	using Upp::vreinterpretq_s32_u32;
	using Upp::vreinterpretq_u64_s64;
	using Upp::vreinterpret_u8_s8;
	using Upp::vreinterpretq_u8_s8;
	using Upp::vreinterpret_u32_s8;
	using Upp::vreinterpretq_u32_s32;
	using Upp::vcombine_f32;
	using Upp::vreinterpret_s8_s32;
	using Upp::vdup_n_s32;
	using Upp::vreinterpret_s32_s8;
	using Upp::vzip_s8;
	using Upp::vcombine_s8;
	using Upp::vreinterpret_u8_u32;
	using Upp::vreinterpret_u32_u8;
	using Upp::vzip_u8;
	using Upp::vcombine_u8;
	using Upp::vreinterpret_u32_s16;
	using Upp::vreinterpret_s16_u32;
	using Upp::vzip_u32;
	using Upp::vzip_s16;
	using Upp::vcombine_s16;
	using Upp::vzip_u16;
	using Upp::vcombine_u16;
	using Upp::vcombine_s32;
	using Upp::vcombine_u32;
	using Upp::vzip_u32;
	using Upp::vcombine_s8;
	using Upp::vcombine_u8;
	using Upp::vcombine_s16;
	using Upp::vcombine_u16;
	using Upp::int8x8x2_t;
	using Upp::uint8x8x2_t;
	using Upp::int16x4x2_t;
	using Upp::uint16x4x2_t;
	using Upp::vreinterpret_u32_u16;
	using Upp::vreinterpret_u16_u32;
	using Upp::vrev64_f32;
	using Upp::vrev64q_f32;
	using Upp::vcombine_f32;
	using Upp::vget_high_f32;
	using Upp::vget_low_f32;
	using Upp::vrev64_s8;
	using Upp::vreinterpret_s8_s32;
	using Upp::vdup_n_s32;
	using Upp::vreinterpret_s32_s8;
	using Upp::vrev64q_s8;
	using Upp::vcombine_s8;
	using Upp::vget_high_s8;
	using Upp::vget_low_s8;
	using Upp::vrev64_u8;
	using Upp::vreinterpret_u8_u32;
	using Upp::vreinterpret_u32_u8;
	using Upp::vrev64q_u8;
	using Upp::vcombine_u8;
	using Upp::vget_high_u8;
	using Upp::vget_low_u8;
	using Upp::vrev64_s16;
	using Upp::vrev64q_s16;
	using Upp::vcombine_s16;
	using Upp::vget_high_s16;
	using Upp::vget_low_s16;
	using Upp::vrev64_u16;
	using Upp::vrev64q_u16;
	using Upp::vcombine_u16;
	using Upp::vget_high_u16;
	using Upp::vget_low_u16;
	using Upp::vrev64_s32;
	using Upp::vrev64q_s32;
	using Upp::vcombine_s32;
	using Upp::vget_high_s32;
	using Upp::vget_low_s32;
	using Upp::vrev64_u32;
	using Upp::vrev64q_u32;
	using Upp::vcombine_u32;
	using Upp::vget_high_u32;
	using Upp::vget_low_u32;
	using Upp::vcombine_s64;
	using Upp::vget_high_s64;
	using Upp::vget_low_s64;
	using Upp::vcombine_u64;
	using Upp::vget_high_u64;
	using Upp::vget_low_u64;
	using Upp::vabs_f32;
	using Upp::vabsq_f32;
	using Upp::vabs_s8;
	using Upp::vabsq_s8;
	using Upp::vabs_s16;
	using Upp::vabsq_s16;
	using Upp::vabs_s32;
	using Upp::vabsq_s32;
	using Upp::vabsq_s64;
	using Upp::vreinterpretq_s32_f32;
	using Upp::vreinterpret_s32_f32;
	using Upp::vreinterpretq_s32_f32;
	using Upp::vreinterpret_s32_f32;
	using Upp::vaddv_f32;
	using Upp::vaddvq_f32;
	using Upp::vpadd_f32;
	using Upp::vadd_f32;
	using Upp::vadd_s8;
	using Upp::vpadd_s8;
	using Upp::vpadd_u8;
	using Upp::vadd_u8;
	using Upp::vadd_s16;
	using Upp::vpadd_s16;
	using Upp::vadd_u16;
	using Upp::vpadd_u16;
	using Upp::vadd_s32;
	using Upp::vpadd_s32;
	using Upp::vadd_u32;
	using Upp::vpadd_u32;
	using Upp::vminv_f32;
	using Upp::vminvq_f32;
	using Upp::vpmin_f32;
	using Upp::vmin_f32;
	using Upp::vpmin_s8;
	using Upp::vmin_s8;
	using Upp::vminv_s8;
	using Upp::vminvq_s8;
	using Upp::vpmin_u8;
	using Upp::vmin_u8;
	using Upp::vminv_u8;
	using Upp::vminvq_u8;
	using Upp::vpmin_s16;
	using Upp::vmin_s16;
	using Upp::vminv_s16;
	using Upp::vminvq_s16;
	using Upp::vpmin_u16;
	using Upp::vmin_u16;
	using Upp::vminv_u16;
	using Upp::vminvq_u16;
	using Upp::vpmin_s32;
	using Upp::vmin_s32;
	using Upp::vminv_s32;
	using Upp::vminvq_s32;
	using Upp::vpmin_u32;
	using Upp::vmin_u32;
	using Upp::vminv_u32;
	using Upp::vminvq_u32;
	using Upp::vmaxv_f32;
	using Upp::vmaxvq_f32;
	using Upp::vpmax_f32;
	using Upp::vmax_f32;
	using Upp::vpmax_s8;
	using Upp::vmax_s8;
	using Upp::vmaxv_s8;
	using Upp::vmaxvq_s8;
	using Upp::vpmax_u8;
	using Upp::vmax_u8;
	using Upp::vmaxv_u8;
	using Upp::vmaxvq_u8;
	using Upp::vpmax_s16;
	using Upp::vmax_s16;
	using Upp::vmaxv_s16;
	using Upp::vmaxvq_s16;
	using Upp::vpmax_u16;
	using Upp::vmax_u16;
	using Upp::vmaxv_u16;
	using Upp::vmaxvq_u16;
	using Upp::vpmax_s32;
	using Upp::vmax_s32;
	using Upp::vmaxv_s32;
	using Upp::vmaxvq_s32;
	using Upp::vpmax_u32;
	using Upp::vmax_u32;
	using Upp::vmaxv_u32;
	using Upp::vmaxvq_u32;
	using Upp::vaddvq_s64;
	using Upp::vaddvq_u64;
	using Upp::vaddv_s8;
	using Upp::vaddvq_s8;
	using Upp::vaddv_s16;
	using Upp::vaddvq_s16;
	using Upp::vaddv_u8;
	using Upp::vaddvq_u8;
	using Upp::vaddv_u16;
	using Upp::vaddvq_u16;
	using Upp::vaddv_s32;
	using Upp::vaddvq_s32;
	using Upp::vaddv_u32;
	using Upp::vaddvq_u32;
	using Upp::vreinterpret_f32_s32;
	using Upp::vreinterpretq_f32_s32;
	using Upp::vmul_f32;
	using Upp::vget_low_f32;
	using Upp::vget_high_f32;
	using Upp::vreinterpret_s8_s32;
	using Upp::vdup_n_s32;
	using Upp::vmul_s8;
	using Upp::vrev16_s8;
	using Upp::vrev32_s8;
	using Upp::vget_low_s8;
	using Upp::vget_high_s8;
	using Upp::vreinterpret_u8_u32;
	using Upp::vdup_n_u32;
	using Upp::vmul_u8;
	using Upp::vrev16_u8;
	using Upp::vrev32_u8;
	using Upp::vget_low_u8;
	using Upp::vget_high_u8;
	using Upp::vmul_s16;
	using Upp::vrev32_s16;
	using Upp::vget_low_s16;
	using Upp::vget_high_s16;
	using Upp::vmul_u16;
	using Upp::vrev32_u16;
	using Upp::float32x2x2_t;
	using Upp::float32x4x2_t;
	using Upp::int8x8x2_t;
	using Upp::int8x16x2_t;
	using Upp::uint8x8x2_t;
	using Upp::uint8x16x2_t;
	using Upp::int32x2x2_t;
	using Upp::int32x4x2_t;
	using Upp::uint32x2x2_t;
	using Upp::uint32x4x2_t;
	using Upp::int16x8x2_t;
	using Upp::int16x4x2_t;
	using Upp::uint16x8x2_t;
	using Upp::vzip_f32;
	using Upp::vzipq_f32;
	using Upp::vzip_s8;
	using Upp::vzipq_s8;
	using Upp::vzip_u8;
	using Upp::vzipq_u8;
	using Upp::vzip_s32;
	using Upp::vzipq_s32;
	using Upp::vzip_u32;
	using Upp::vzipq_u32;
	using Upp::vzip_s16;
	using Upp::vzipq_s16;
	using Upp::vzip_u16;
	using Upp::vzipq_u16;
	using Upp::vreinterpret_s16_s8;
	using Upp::vreinterpret_s32_s16;
	using Upp::vreinterpret_u16_u8;
	using Upp::vreinterpret_u32_u8;
	using Upp::vzip1q_s64;
	using Upp::vzip2q_s64;
	using Upp::vzip1q_u64;
	using Upp::vzip2q_u64;
	using Upp::vget_low_s64;
	using Upp::vget_high_s64;
	using Upp::vget_low_u64;
	using Upp::vget_high_u64;
	using Upp::vcombine_s64;
	using Upp::vcombine_u64;
	using Upp::vbsl_f32;
	using Upp::vbslq_f32;
	using Upp::vbsl_s8;
	using Upp::vbslq_s8;
	using Upp::vbsl_u8;
	using Upp::vbslq_u8;
	using Upp::vbsl_s16;
	using Upp::vbslq_s16;
	using Upp::vbsl_u16;
	using Upp::vbslq_u16;
	using Upp::vbsl_s32;
	using Upp::vbslq_s32;
	using Upp::vbsl_u32;
	using Upp::vbslq_u32;
	using Upp::vbslq_s64;
	using Upp::vbslq_u64;
	using Upp::vreinterpret_u32_f32;
	using Upp::vreinterpretq_u32_f32;
	using Upp::vreinterpret_u8_s8;
	using Upp::vreinterpretq_u8_s8;
	using Upp::vreinterpret_u16_s16;
	using Upp::vreinterpretq_u16_s16;
	using Upp::vreinterpret_u32_s32;
	using Upp::vreinterpretq_u32_s32;
	using Upp::vreinterpretq_u64_s64;
	using Upp::vreinterpret_u8_u32;
	using Upp::vreinterpret_u32_u8;
	using Upp::vrndn_f32;
	using Upp::vrndnq_f32;
	using Upp::vrndm_f32;
	using Upp::vrndmq_f32;
	using Upp::vrndp_f32;
	using Upp::vrndpq_f32;
	using Upp::vrnda_f32;
	using Upp::vrndaq_f32;
	using Upp::vrnd_f32;
	using Upp::vrndq_f32;
	using Upp::vdup_n_u8;
	using Upp::vdupq_n_u8;
	using Upp::vorr_u8;
	using Upp::vorrq_u8;
	using Upp::vcge_u8;
	using Upp::vcgeq_u8;
	using Upp::vmul_u8;
	using Upp::vdup_n_u16;
	using Upp::vdupq_n_u16;
	using Upp::vorr_u16;
	using Upp::vorrq_u16;
	using Upp::vcge_u16;
	using Upp::vcgeq_u16;
	using Upp::vmul_u16;
	using Upp::vdup_n_u32;
	using Upp::vdupq_n_u32;
	using Upp::vorr_u32;
	using Upp::vorrq_u32;
	using Upp::vcge_u32;
	using Upp::vcgeq_u32;
	using Upp::vmul_u32;
	using Upp::vrsqrteq_f32;
	using Upp::vrsqrtsq_f32;
	using Upp::vmulq_f32;
	using Upp::vrsqrte_f32;
	using Upp::vrsqrts_f32;
	using Upp::vmul_f32;
	using Upp::vrecpeq_f32;
	using Upp::vrecpsq_f32;
	using Upp::vrecpe_f32;
	using Upp::vrecps_f32;
	using Upp::vsqrtq_f32;
	using Upp::vsqrt_f32;
	using Upp::vdivq_f32;
	using Upp::vdiv_f32;
	using Upp::vmovn_u32;
	using Upp::vmovl_u16;
	using Upp::vaddq_u8;
	using Upp::vaddq_u16;
	using Upp::vaddq_u32;
	using Upp::vandq_u8;
	using Upp::vandq_u16;
	using Upp::vandq_u32;
	using Upp::vbicq_u8;
	using Upp::vclzq_u8;
	using Upp::vclzq_u16;
	using Upp::vclzq_u32;
	using Upp::vcntq_u8;
	using Upp::veorq_u8;
	using Upp::veorq_u16;
	using Upp::veorq_u32;
	using Upp::vmaxq_u8;
	using Upp::vmaxq_u16;
	using Upp::vmaxq_u32;
	using Upp::vminq_u8;
	using Upp::vminq_u16;
	using Upp::vminq_u32;
	using Upp::vmovl_u8;
	using Upp::vmovl_u16;
	using Upp::vmovl_u32;
	using Upp::vmovn_u16;
	using Upp::vmovn_u32;
	using Upp::vmovn_u64;
	using Upp::vmvnq_u8;
	using Upp::vmvnq_u16;
	using Upp::vmvnq_u32;
	using Upp::vorrq_u8;
	using Upp::vorrq_u16;
	using Upp::vorrq_u32;
	using Upp::vqaddq_u8;
	using Upp::vqaddq_u16;
	using Upp::vqaddq_u32;
	using Upp::vqmovn_u16;
	using Upp::vqmovn_u32;
	using Upp::vqmovn_u64;
	using Upp::vqsubq_u8;
	using Upp::vqsubq_u16;
	using Upp::vqsubq_u32;
	using Upp::vrecpeq_f32;
	using Upp::vrecpsq_f32;
	using Upp::vrev16q_u8;
	using Upp::vrev32q_u16;
	using Upp::vrev32q_u8;
	using Upp::vrev64q_u16;
	using Upp::vrev64q_u32;
	using Upp::vrev64q_u8;
	using Upp::vrhaddq_u8;
	using Upp::vrhaddq_u16;
	using Upp::vrhaddq_u32;
	using Upp::vshlq_u8;
	using Upp::vshlq_u16;
	using Upp::vshlq_u32;
	using Upp::vshlq_u64;
	using Upp::vsubq_u8;
	using Upp::vsubq_u16;
	using Upp::vsubq_u32;
	using Upp::vdupq_n_f64;
	using Upp::vaddq_f64;
	using Upp::vsubq_f64;
	using Upp::vdupq_n_f64;
	using Upp::vaddq_f64;
	using Upp::vsubq_f64;
	using Upp::vnegq_f64;
	using Upp::vmulq_f64;
	using Upp::vdivq_f64;
	using Upp::vfmaq_f64;
	using Upp::vfmsq_f64;
	using Upp::vmlaq_f64;
	using Upp::vmlsq_f64;
	using Upp::vminq_f64;
	using Upp::vminnmq_f64;
	using Upp::vmaxnmq_f64;
	using Upp::vmaxq_f64;
	using Upp::vreinterpretq_f64_u64;
	using Upp::vandq_u64;
	using Upp::vreinterpretq_u64_f64;
	using Upp::vorrq_u64;
	using Upp::veorq_u64;
	using Upp::vbicq_u64;
	using Upp::vcleq_f64;
	using Upp::vcltq_f64;
	using Upp::vmvnq_u32;
	using Upp::vreinterpretq_f64_u32;
	using Upp::vreinterpretq_u32_u64;
	using Upp::vcgeq_f64;
	using Upp::vceqq_f64;
	using Upp::vcombine_f64;
	using Upp::vget_high_f64;
	using Upp::vget_low_f64;
	using Upp::vabsq_f64;
	using Upp::vreinterpretq_f64_s64;
	using Upp::vaddvq_f64;
	using Upp::vmul_f64;
	using Upp::vminvq_f64;
	using Upp::vmaxvq_f64;
	using Upp::vzip1q_f64;
	using Upp::vzip2q_f64;
	using Upp::vbslq_f64;
	using Upp::vreinterpretq_u64_f64;
	using Upp::vrndnq_f64;
	using Upp::vrndmq_f64;
	using Upp::vrndpq_f64;
	using Upp::vrndaq_f64;
	using Upp::vrndq_f64;
	using Upp::vdupq_n_u64;
	using Upp::vrsqrteq_f64;
	using Upp::vsqrtq_f64;
	using Upp::vreinterpretq_s64_f64;
	using Upp::vreinterpretq_f64_s64;
	using Upp::float32_t;
	using Upp::float64_t;
	using Upp::vcvtq_s64_f64;
	using Upp::vcvt_f64_f32;
	using Upp::vcvtq_u64_f64;
	using Upp::vmovl_s32;
	using Upp::vget_low_s32;
	using Upp::vcvtq_s32_f32;
	using Upp::vcvt_s32_f32;
	using Upp::vmovl_u32;
	using Upp::vget_low_u32;
	using Upp::vcvtq_u32_f32;
	using Upp::vcvt_u32_f32;
	using Upp::vmovn_s32;
	using Upp::vcombine_s16;
	using Upp::vmovn_u32;
	using Upp::vcombine_u16;
	using Upp::vmovn_u32;
	using Upp::vcombine_u32;
	using Upp::vcombine_s32;
	using Upp::vmovn_s16;
	using Upp::vcombine_s8;
	using Upp::vcombine_s16;
	using Upp::vmovn_s16;
	using Upp::vreinterpret_s32_s8;
	using Upp::vcvtq_f32_s32;
	using Upp::vmovl_s16;
	using Upp::vget_low_s16;
	using Upp::vmovl_s8;
	using Upp::vget_low_s8;
	using Upp::vreinterpret_s8_s32;
	using Upp::vdup_n_s32;
	using Upp::vcvt_f32_s32;
	using Upp::vget_low_s32;
	using Upp::vmovl_s32;
	using Upp::vreinterpret_s8_s32;
	using Upp::vget_low_u16;
	using Upp::vmovl_u16;
	using Upp::vmovl_u8;
	using Upp::vget_low_u8;
	using Upp::vreinterpret_u8_u32;
	using Upp::vdup_n_u32;
	using Upp::vget_low_u32;
	using Upp::vmovl_u32;
	using Upp::vget_low_s8;
	using Upp::vmovl_s8;
	using Upp::vget_low_s16;
	using Upp::vmovl_s16;
	using Upp::vget_low_s32;
	using Upp::vmovl_s32;
	using Upp::vget_low_u8;
	using Upp::vmovl_u8;
	using Upp::vget_low_u16;
	using Upp::vmovl_u16;
	using Upp::vget_low_u32;
	using Upp::vmovl_u32;
	using Upp::vget_low_s8;
	using Upp::vmovl_s8;
	using Upp::vcvt_f32_u32;
	using Upp::vcvtq_f32_u32;
	using Upp::vcvtq_f64_s64;
	using Upp::vcvt_f64_s64;
	using Upp::vcvt_f32_f64;
	using Upp::vmovn_s64;
	using Upp::vcvtq_f64_u64;
	using Upp::vreinterpretq_f64_s32;
	using Upp::vreinterpretq_s32_f64;
	using Upp::vreinterpretq_u32_f64;
	using Upp::vreinterpretq_f32_u64;
	using Upp::vfmaq_n_f32;
	using Upp::vfmaq_n_f64;
}
#endif 

#undef Success  
#include <plugin/eigen/Eigen/Dense>
#include <plugin/eigen/unsupported/Eigen/NonLinearOptimization>
#undef Complex
#include <plugin/eigen/unsupported/Eigen/FFT>
#include <plugin/eigen/unsupported/Eigen/CXX11/Tensor>

#include <Functions4U/Defs.h>
#include "MultiDimMatrix.h"

namespace Upp {

template <class T>
using UVector = Upp::Vector<T>;

template <class T>
using UArray = Upp::Array<T>;

template <class T>
using UIndex = Upp::Index<T>;

template<typename _Scalar, ptrdiff_t nx = Eigen::Dynamic, ptrdiff_t ny = Eigen::Dynamic>
struct NonLinearOptimizationFunctor {
	typedef _Scalar Scalar;
	enum {
		InputsAtCompileTime = nx,
		ValuesAtCompileTime = ny
	};
	typedef Eigen::Matrix<double, InputsAtCompileTime, 1> InputType;
	typedef Eigen::Matrix<double, ValuesAtCompileTime, 1> ValueType;
	typedef Eigen::Matrix<double, ValuesAtCompileTime, InputsAtCompileTime> JacobianType;
	
	Eigen::Index unknowns, datasetLen;
	
	NonLinearOptimizationFunctor() : unknowns(InputsAtCompileTime), datasetLen(ValuesAtCompileTime) {}
	NonLinearOptimizationFunctor(int unknowns, int datasetLen) : unknowns(unknowns), datasetLen(datasetLen) {}
	
	ptrdiff_t inputs() const {return ptrdiff_t(unknowns);}
	ptrdiff_t values() const {return ptrdiff_t(datasetLen);}
	virtual void operator() (const InputType& , ValueType* , JacobianType*  = 0) const {};
};

struct BasicNonLinearOptimizationFunctor : NonLinearOptimizationFunctor<double> {
	BasicNonLinearOptimizationFunctor(Function <int(const Eigen::VectorXd &b, Eigen::VectorXd &err)> _function) : function(_function) {}
	int operator()(const Eigen::VectorXd &b, Eigen::VectorXd &fvec) const {return function(b, fvec);}
	Function <int(const Eigen::VectorXd &b, Eigen::VectorXd &err)> function;
};

bool NonLinearOptimization(Eigen::VectorXd &y, Eigen::Index numData, 
			Function <int(const Eigen::VectorXd &y, Eigen::VectorXd &residual)>residual,
			double xtol = Null, double ftol = Null, int maxfev = Null);
bool NonLinearOptimization(Eigen::VectorXd &y, Eigen::Index numData, 
			Function <int(const Eigen::VectorXd &b, Eigen::VectorXd &residual)> Residual,
			double xtol, double ftol, int maxfev, int &ret);
String NonLinearOptimizationError(int error);

bool SolveNonLinearEquations(Eigen::VectorXd &y, Function <int(const Eigen::VectorXd &b, Eigen::VectorXd &residual)> Residual,
			double xtol = Null, int maxfev = Null, double factor = Null);
double SolveNonLinearEquation(double y, Function <double(double b)> Residual, double xtol = Null, int maxfev = Null, double factor = Null);

Eigen::Matrix3d SkewSymmetricMatrix(const Eigen::Vector3d& r);		// Antisymmetric matrix from 3d vector
	
template <class T>
void Xmlize(XmlIO &xml, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &mat) {
	Size_<int64> sz(mat.cols(), mat.rows());
	xml ("size", sz);
	if(xml.IsStoring()) {
		for(int r = 0; r < mat.rows(); r++)
			for(int c = 0; c < mat.cols(); c++) {
				XmlIO io = xml.Add("item");
				T data = mat(r, c);
				Xmlize(io, data);
			}
	} else {
		mat.resize(ptrdiff_t(sz.cy), ptrdiff_t(sz.cx));
		int r = 0, c = 0;
		for(int i = 0; i < xml->GetCount(); i++) 
			if(xml->Node(i).IsTag("item")) {
				XmlIO io = xml.At(i);
				T data;
				Xmlize(io, data);
				mat(r, c) = data;
				++c;
				if (c == sz.cx) {
					c = 0;
					r++;
				}
			}
	}
}

template <class T>
void Xmlize(XmlIO &xml, Eigen::Matrix<T, Eigen::Dynamic, 1> &vec) {
	int64 sz = vec.size();
	xml ("size", sz);
	if(xml.IsStoring()) {
		for(int r = 0; r < sz; r++) {
			XmlIO io = xml.Add("item");
			T data = vec(r);
			Xmlize(io, data);
		}
	} else {
		vec.resize(ptrdiff_t(sz));
		int r = 0;
		for(int i = 0; i < xml->GetCount(); i++)
			if(xml->Node(i).IsTag("item")) {
				XmlIO io = xml.At(i);
				T data;
				Xmlize(io, data);
				vec(r++) = data;
			}
	}
}

template <typename T, int NumIndices>
void Jsonize(JsonIO &io, Eigen::Tensor<T, NumIndices> &mat) {
	Array<T> vector;
	Vector<int64> vsz(NumIndices);
	for (int i = 0; i < NumIndices; ++i)
		vsz[i] = mat.dimension(i);
	io("size", vsz);
	if(io.IsStoring()) {
		vector.SetCount(int(mat.size()));
		Copy(mat, vector);
		io("data", vector);
	} else {
		io("data", vector);
		Eigen::array<Eigen::Index, NumIndices> dims;
		for (int i = 0; i < NumIndices; ++i) 	
			dims[i] = static_cast<Eigen::Index>(vsz[i]);
		mat.resize(dims);
		std::copy(vector.begin(), vector.end(), mat.data());
	}
}

template <class T>
void Jsonize(JsonIO &io, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &mat) {
	Array<T> vector;
	Size_<int64> sz(mat.cols(), mat.rows());
	io("size", sz);
	if(io.IsStoring()) {
		vector.SetCount(int(sz.cx)*int(sz.cy));
		Copy(mat, vector);
		io("data", vector);
	} else {
		io("data", vector);
		mat.resize(ptrdiff_t(sz.cy), ptrdiff_t(sz.cx));
		std::copy(vector.begin(), vector.end(), mat.data());
	}
}

template <class T>
void Jsonize(JsonIO &io, Eigen::Matrix<T, Eigen::Dynamic, 1> &vec) {
	Array<T> vector;
	int64 sz = vec.size();
	io("size", sz);
	if(io.IsStoring()) {
		vector.SetCount(int(sz));
		Copy(vec, vector);
		io("data", vector);
	} else {
		io("data", vector);
		vec.resize(ptrdiff_t(sz));
		Copy(vector, vec);
	}
}

template <class T>
void Serialize(Stream& stream, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &mat) {
	Size_<int64> sz(mat.cols(), mat.rows());
	stream % sz;
	if(stream.IsStoring()) {
		for(int r = 0; r < mat.rows(); r++)
			for(int c = 0; c < mat.cols(); c++) {
				T data = mat(r, c);
				stream % data;
			}
	} else {
		mat.resize(ptrdiff_t(sz.cy), ptrdiff_t(sz.cx));
		int r = 0, c = 0;
		for(int i = 0; i < sz.cy*sz.cx; i++) {
			T data;
			stream % data;
			mat(r, c) = data;
			++c;
			if (c == sz.cx) {
				c = 0;
				r++;
			}
			if (r == sz.cy)
				break;
		}
	}
}

template <class T>
void Serialize(Stream& stream, Eigen::Matrix<T, Eigen::Dynamic, 1> &vec) {
	int64 sz = vec.size();
	stream % sz;
	if(stream.IsStoring()) {
		for (int i = 0; i < sz; ++i) {
			T data = vec(i);
			stream % data;
		}
	} else {
		vec.resize(ptrdiff_t(sz));
		for (int i = 0; i < sz; ++i) {
			T data;
			stream % data;
			vec(i) = data;
		}
	}
}

// These functions serve both for Eigen, std and U++ Vectors

template <class Range>
void Resize(Range &v, size_t len) {v.SetCount(int(len));}
template <class Range>
void Resize(Range &v, size_t len, const typename Range::value_type& init) {
	v.SetCount(int(len));
	std::fill(v.begin(), v.end(), init);
}

template <class Range>
void ResizeConservative(Range &v, size_t len) {v.SetCount(int(len));}
template <class Range>
void ResizeConservative(Range &v, size_t len, const typename Range::value_type& init) {v.SetCount(int(len), init);}
template <class Range>
void Clear(Range &v) {v.Clear();}

template <typename T>
void Resize(Eigen::Matrix<T, Eigen::Dynamic, 1> &v, size_t len) {v.resize((Eigen::Index)len);}
template <typename T>
void Resize(Eigen::Matrix<T, Eigen::Dynamic, 1> &v, size_t len, const T& init) {v.setConstant(len, 1, init);}
template <typename T>
void ResizeConservative(Eigen::Matrix<T, Eigen::Dynamic, 1> &v, size_t len) {v.conservativeResize(len);}
template <typename T>
void ResizeConservative(Eigen::Matrix<T, Eigen::Dynamic, 1> &v, size_t len, const T& init) {
	size_t len0 = v.size();
	v.conservativeResize(len);
	if (len > len0)
		std::fill(&v[len0], v.data() + len, init);
}
template <typename T>
void Clear(Eigen::Matrix<T, Eigen::Dynamic, 1> &v) 				{v.resize(0);}
template <typename T>
void Clear(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &m) {m.resize(0, 0);}
template <typename T, int NumIndices_>
void Clear(Eigen::Tensor<T, NumIndices_> &m) 					{m = Eigen::Tensor<T, NumIndices_>();}

template <typename T>
void PrePad(Eigen::Matrix<T, Eigen::Dynamic, 1> &v, size_t len, const T& init) {
	size_t len0 = v.size();
	v.conservativeResize(len);
	if (len > len0) {
		size_t delta = len - len0;
		std::copy(&v[len0 - delta], v.data() + len0, &v[len0]);
		std::copy(v.data(), v.data() + len0 - delta, &v[delta]);
		std::fill(v.data(), v.data() + delta, init);
	}
}


template <typename T>
void Resize(std::vector<T> &v, size_t len) {v.resize(len);}
template <typename T>
void Resize(std::vector<T> &v, size_t len, const T& init) {
	v.resize(len);
	std::fill(v.begin(), v.end(), init);
}
template <typename T>
void ResizeConservative(std::vector<T> &v, size_t len) {v.resize(len);}
template <typename T>
void ResizeConservative(std::vector<T> &v, size_t len, const T& init) {v.resize(len, init);}
template <typename T>
void Clear(std::vector<T> &v) {v.clear();}

#define PostPad ResizeConservative


template <typename T>
auto Begin(const std::vector<T> &v)		{return v.begin();}
template <typename T>
auto Begin(std::vector<T> &v)			{return v.begin();}
template <typename T>
auto End(const std::vector<T> &v)		{return v.end();}
template <typename T>
auto End(std::vector<T> &v)				{return v.end();}

template <typename T>
auto Begin(const Eigen::Matrix<T, Eigen::Dynamic, 1> &v){return v.data();}
template <typename T>
auto Begin(Eigen::Matrix<T, Eigen::Dynamic, 1> &v)		{return v.data();}
template <typename T>
auto End(const Eigen::Matrix<T, Eigen::Dynamic, 1> &v)	{return v.data() + v.size();}

template <typename T>
auto Begin(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &v)	{return v.data();}
template <typename T>
auto Begin(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &v)			{return v.data();}
template <typename T>
auto End(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &v)		{return v.data() + v.size();}

template <typename T, int NumIndices>
auto Begin(const Eigen::Tensor<T, NumIndices> &v)		{return v.data();}
template <typename T, int NumIndices>
auto Begin(Eigen::Tensor<T, NumIndices> &v)				{return v.data();}
template <typename T, int NumIndices>
auto End(const Eigen::Tensor<T, NumIndices> &v)			{return v.data() + v.size();}

template <class Range>
auto Begin(const Range &v)				{return v.Begin();}
template <class Range>
auto Begin(Range &v)					{return v.Begin();}
template <class Range>
auto End(const Range &v)				{return v.End();}
template <class Range>
auto End(Range &v)						{return v.End();}


template <class Range>
auto &First(Range &data) {return data[0];}

template <class Range>
auto &Last(Range &data) {return data[data.size()-1];}

template <typename T>			// To avoid problem with Upp
void ReverseX(std::vector<T> &v) {std::reverse(v.begin(), v.end());}

template <typename T>			// To avoid problem with Upp
void ReverseX(Eigen::Matrix<T, Eigen::Dynamic, 1> &v) {v.reverseInPlace();}

template <typename T>			// To avoid problem with Upp
void ReverseX(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &v) {v.reverseInPlace();}

template <class Range>
void ReverseX(Range &v) {		// To avoid problem with Upp
	typename Range::value_type *first = v.begin();
	typename Range::value_type *last = v.end();
	while ((first != last) && (first != --last)) 
		Swap(*first++, *last);
}

template <typename T>			
void Rotate(std::vector<T> &v, int shift) {
	if (shift > 0)
		std::rotate(v.begin(), v.begin() + shift, v.end());
	else if (shift < 0)
		std::rotate(v.rbegin(), v.rbegin() - shift, v.rend());
}

template <class Range>			
void Rotate(Range &v, int k) {
	auto rotate = [](Range &v, int start, int end) {
	    while (start < end) {
	        Swap(v[start], v[end]);
	        start++;
	        end--;
	    }
	};
	int n = v.size();
	k = k % n;  // Handle cases where k is greater than the array size
	if (k > 0) {
	    rotate(v, 0, n - 1);
	    rotate(v, 0, k - 1);
	    rotate(v, k, n - 1);
	} else if (k > 0) {
	    rotate(v, 0, k - 1);
	    rotate(v, k, n - 1);
	    rotate(v, 0, n - 1);
	}
}

template <class Range>
void CopyRowMajor(Range &in, int nrows, int ncols, Eigen::Matrix<typename Range::value_type, Eigen::Dynamic, Eigen::Dynamic> &out) {
	out = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(Begin(in), nrows, ncols);
}

template <typename T>
void CopyRowMajor(T *in, int nrows, int ncols, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &out) {
	out = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(in, nrows, ncols);
}

template <class Range>
void CopyRowMajor(const Eigen::Matrix<typename Range::value_type, Eigen::Dynamic, Eigen::Dynamic> &in, Range &out) {
	Resize(out, in.size());
	typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMajMat;
	RowMajMat::Map(Begin(out), in.rows(), in.cols()) = in;
}

template <class Range1, class Range2>
void Copy(const Range1& in, Range2 &out) {
	Resize(out, (size_t)in.size());
	std::copy(Begin(in), End(in), Begin(out));
}

template <class Range1, class Range2>
void AppendX(const Range1& in, Range2 &out) {		// Append is already used in Core/Obsolete
	size_t outsize = (size_t)out.size();
	ResizeConservative(out, (size_t)in.size() + outsize);
	std::copy(Begin(in), End(in), Begin(out) + outsize);
}

template <class Range1, class Range2>
void Block(Range1& in, Range2& out, int begin, int len) {
	ASSERT(begin < in.size() && len > 0 && begin + len <= in.size());
	Resize(out, len);
	std::copy(in + begin, in + begin + len, Begin(out));
}


template <class T>
void Swap(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &A, int rc1, int rc2) {
 	A.row(rc1).swap(A.row(rc2));	// Swap rows rc1 and rc2
 	A.col(rc1).swap(A.col(rc2));	// Swap columns rc1 and rc2
}

template <class T>
void Swap(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &A1, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &A2, int rc1, int rc2) {
	ASSERT(A1.rows() == A2.rows() && A1.cols() == A2.cols());
	
	A1.row(rc1).swap(A2.row(rc2));
    A1.col(rc1).swap(A2.col(rc2));	
}

template <typename T>
inline Eigen::Matrix<T, Eigen::Dynamic, 1> Segment(const Eigen::Matrix<T, Eigen::Dynamic, 1> &d, int ifrom, int num = -1) {
	if (num < 0)
		num = (int)d.size()- ifrom;
	return d.segment(ifrom, num);
}

template <typename T>
inline std::vector<T> Segment(const std::vector<T> &d, int ifrom, int num = -1) {
	if (num < 0)
		num = d.size()- ifrom;
	return std::vector<T>(d.begin() + ifrom, d.begin() + ifrom + num);
}

template <class Range>
inline Range Segment(const Range &d, int ifrom, int num = -1) {
	Range a;
	if (num < 0)
		num = d.size()- ifrom;
	else {
		if (ifrom + num >= d.size()) {
			num = d.size() - ifrom;
			if (num <= 0)
				return a;
		}
	}
	Resize(a, num);
	std::copy(Begin(d) + ifrom, Begin(d) + ifrom + num, Begin(a));
	return a;
}

template <class Range>
inline Range Mid(const Range &d, int ifrom, int num) {return Segment(d, ifrom, num);}

template <class Range>
inline Range Left(const Range &d, int num) {
	Range a;
	if (num >= d.size()) {
		a = clone(d);
		return a; 
	}
	Resize(a, num);
	std::copy(Begin(d), Begin(d) + num, Begin(a));
	return a;
}

template <class Range>
inline Range Right(const Range &d, int num) {
	Range a;
	if (num >= d.size()) {
		a = clone(d);
		return a; 
	}
	Resize(a, num);
	std::copy(Begin(d) + d.size() - num, Begin(d) + d.size(), Begin(a));
	return a;
}

template <typename T>
inline void Remove(Vector<T> &d, int id) {d.Remove(id);}

template <typename T>
inline void Remove(Eigen::Matrix<T, Eigen::Dynamic, 1> &d, int id) {
	Eigen::Index sz = d.size();
	Eigen::Index right = sz-id-1;
	d.segment(id, right) = d.segment(id+1, right);
	d.conservativeResize(sz-1);
}

template <class Range>
inline void Remove(Range &d, int id) {
	d.erase(d.begin() + id);
}

template <class T>
bool IsNum(const Eigen::Matrix<T, Eigen::Dynamic, 1> &a) {
	if (a.size() == 0)
		return false;
	for (int i = 0; i < a.size(); ++i) {
		if (!IsNum(a(i)))
			return false;
	}
	return true;
}

template <class T>
bool IsNum(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &a) {
	if (a.size() == 0)
		return false;
	for (int i = 0; i < a.size(); ++i) {
		if (!IsNum(a.array()(i)))
			return false;
	}
	return true;
}

template <class T>
void Nvl2(Eigen::Matrix<T, Eigen::Dynamic, 1> &a, T val) {
	if (a.size() == 0)
		return;
	for (int i = 0; i < a.size(); ++i)
		if (!IsNum(a(i)))
			a(i) = val;
}

template <class T>
void Nvl2(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &a, T val) {
	if (a.size() == 0)
		return;
	for (int i = 0; i < a.size(); ++i)
		if (!IsNum(a.array()(i)))
			a(i) = val;
}

#define EigenNull	Eigen::MatrixXd()

template<typename T>
using  MatrixType = Eigen::Matrix<T,Eigen::Dynamic, Eigen::Dynamic>;
   
template<typename Scalar, int rank, typename sizeType>
auto TensorToMatrix(const Eigen::Tensor<Scalar,rank> &tensor, const sizeType rows, const sizeType cols) { 
    return Eigen::Map<const MatrixType<Scalar>> (tensor.data(), rows, cols);
}

template<typename Scalar, typename... Dims>
auto MatrixToTensor(const MatrixType<Scalar> &matrix, Dims... dims) {
    constexpr int rank = sizeof... (Dims);
    return Eigen::TensorMap<Eigen::Tensor<const Scalar, rank>>(matrix.data(), {dims...});
}

template<typename T>
decltype(auto) TensorLayoutSwap(T&& t) {
	return Eigen::TensorLayoutSwapOp<typename std::remove_reference<T>::type>(t);
}
 
}


namespace Eigen {

template<> struct NumTraits<Upp::Complex> : GenericNumTraits<Upp::Complex>
{
  typedef double Real;
  typedef typename NumTraits<double>::Literal Literal;
  enum {
    IsComplex = 1,
    RequireInitialization = NumTraits<Real>::RequireInitialization,
    ReadCost = 2 * NumTraits<Real>::ReadCost,
    AddCost = 2 * NumTraits<Real>::AddCost,
    MulCost = 4 * NumTraits<Real>::MulCost + 2 * NumTraits<Real>::AddCost
  };

  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
  static inline Real epsilon() { return NumTraits<Real>::epsilon(); }
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
  static inline Real dummy_precision() { return NumTraits<Real>::dummy_precision(); }
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
  static inline int digits10() { return NumTraits<Real>::digits10(); }
};

}
#endif