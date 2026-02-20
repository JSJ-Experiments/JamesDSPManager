#include <math.h>
#include <string.h>
#include <float.h>
#include "../jdsp_header.h"

static inline float clampf(float x, float lo, float hi)
{
	if (x < lo) return lo;
	if (x > hi) return hi;
	return x;
}

static void spectrumBiquadClear(SpectrumExtBiquad *bq)
{
	bq->z1[0] = bq->z1[1] = 0.0f;
	bq->z2[0] = bq->z2[1] = 0.0f;
}

static inline float spectrumBiquadProcess(SpectrumExtBiquad *bq, float in, int ch)
{
	float y = bq->b0 * in + bq->z1[ch];
	bq->z1[ch] = bq->b1 * in - bq->a1 * y + bq->z2[ch];
	bq->z2[ch] = bq->b2 * in - bq->a2 * y;
	if (fabsf(bq->z1[ch]) < 1.0e-20f) bq->z1[ch] = 0.0f;
	if (fabsf(bq->z2[ch]) < 1.0e-20f) bq->z2[ch] = 0.0f;
	return y;
}

static void spectrumBiquadSet(SpectrumExtBiquad *bq, float fs, float fc, float q, int highpass)
{
	fc = clampf(fc, 20.0f, fs * 0.49f);
	q = clampf(q, 0.1f, 10.0f);
	float w0 = 2.0f * (float)M_PI * fc / fs;
	float c = cosf(w0);
	float s = sinf(w0);
	float alpha = s / (2.0f * q);
	float b0, b1, b2;
	if (highpass)
	{
		b0 = (1.0f + c) * 0.5f;
		b1 = -(1.0f + c);
		b2 = (1.0f + c) * 0.5f;
	}
	else
	{
		b0 = (1.0f - c) * 0.5f;
		b1 = 1.0f - c;
		b2 = (1.0f - c) * 0.5f;
	}
	float a0 = 1.0f + alpha;
	float a1 = -2.0f * c;
	float a2 = 1.0f - alpha;
	bq->b0 = b0 / a0;
	bq->b1 = b1 / a0;
	bq->b2 = b2 / a0;
	bq->a1 = a1 / a0;
	bq->a2 = a2 / a0;
}

static inline float spectrumHarmonics(float x, const double harmonics[10])
{
	float s = tanhf(x);
	float p = s;
	float out = 0.0f;
	for (int i = 0; i < 10; i++)
	{
		p *= s;
		out += (float)harmonics[i] * p;
	}
	return out;
}

static void SpectrumExtensionRefreshLocked(JamesDSPLib *jdsp)
{
	SpectrumExtension *ext = &jdsp->spectrumExt;
	float fs = jdsp->fs;
	float hpFc = (float)ext->referenceFreq;
	float lpFc = fs * 0.5f - (float)ext->lpOffsetHz;
	if (lpFc < 200.0f)
		lpFc = 200.0f;
	spectrumBiquadSet(&ext->hp, fs, hpFc, ext->hpQ, 1);
	spectrumBiquadSet(&ext->lp, fs, lpFc, ext->lpQ, 0);
	spectrumBiquadClear(&ext->hp);
	spectrumBiquadClear(&ext->lp);
}

void SpectrumExtensionConstructor(JamesDSPLib *jdsp)
{
	SpectrumExtension *ext = &jdsp->spectrumExt;
	memset(ext, 0, sizeof(*ext));
	ext->enabled = 0;
	ext->safetyEnabled = 1;
	ext->strengthLinear = 1.0f;
	ext->referenceFreq = 7600;
	ext->wetMix = 1.0f;
	ext->dryMix = 1.0f;
	ext->postGain = 1.0f;
	ext->hpQ = 0.717f;
	ext->lpQ = 0.717f;
	ext->lpOffsetHz = 2000;
	for (int i = 0; i < 10; i++)
		ext->harmonics[i] = (i % 2 == 0) ? 0.02 : 0.0;
	SpectrumExtensionRefresh(jdsp);
}

void SpectrumExtensionRefresh(JamesDSPLib *jdsp)
{
	jdsp_lock(jdsp);
	SpectrumExtensionRefreshLocked(jdsp);
	jdsp_unlock(jdsp);
}

void SpectrumExtensionSetParam(JamesDSPLib *jdsp, float strengthLinear, int referenceFreq, float wetMix, float postGainDb, char safetyEnabled, float hpQ, float lpQ, int lpOffsetHz, const double harmonics[10])
{
	jdsp_lock(jdsp);
	SpectrumExtension *ext = &jdsp->spectrumExt;
	ext->strengthLinear = clampf(strengthLinear, 0.0f, 3.9810717f);
	ext->referenceFreq = (int)clampf((float)referenceFreq, 800.0f, 20000.0f);
	ext->wetMix = clampf(wetMix, 0.0f, 1.0f);
	ext->dryMix = 1.0f;
	if (postGainDb < -24.0f) postGainDb = -24.0f;
	if (postGainDb > 24.0f) postGainDb = 24.0f;
	ext->postGain = db2magf(postGainDb);
	ext->safetyEnabled = safetyEnabled ? 1 : 0;
	ext->hpQ = clampf(hpQ, 0.1f, 10.0f);
	ext->lpQ = clampf(lpQ, 0.1f, 10.0f);
	ext->lpOffsetHz = (int)clampf((float)lpOffsetHz, 200.0f, 12000.0f);
	for (int i = 0; i < 10; i++)
	{
		double h = harmonics[i];
		if (!isfinite(h))
			h = 0.0;
		if (h < -2.0) h = -2.0;
		if (h > 2.0) h = 2.0;
		ext->harmonics[i] = h;
	}
	SpectrumExtensionRefreshLocked(jdsp);
	jdsp_unlock(jdsp);
}

void SpectrumExtensionEnable(JamesDSPLib *jdsp)
{
	jdsp_lock(jdsp);
	jdsp->spectrumExt.enabled = 1;
	spectrumBiquadClear(&jdsp->spectrumExt.hp);
	spectrumBiquadClear(&jdsp->spectrumExt.lp);
	jdsp_unlock(jdsp);
}

void SpectrumExtensionDisable(JamesDSPLib *jdsp)
{
	jdsp_lock(jdsp);
	jdsp->spectrumExt.enabled = 0;
	spectrumBiquadClear(&jdsp->spectrumExt.hp);
	spectrumBiquadClear(&jdsp->spectrumExt.lp);
	jdsp_unlock(jdsp);
}

void SpectrumExtensionProcess(JamesDSPLib *jdsp, size_t n)
{
	SpectrumExtension *ext = &jdsp->spectrumExt;
	for (size_t i = 0; i < n; i++)
	{
		float inL = jdsp->tmpBuffer[0][i];
		float inR = jdsp->tmpBuffer[1][i];
		float hpL = spectrumBiquadProcess(&ext->hp, inL, 0);
		float hpR = spectrumBiquadProcess(&ext->hp, inR, 1);
		float harmL = spectrumHarmonics(hpL, ext->harmonics) * ext->strengthLinear;
		float harmR = spectrumHarmonics(hpR, ext->harmonics) * ext->strengthLinear;
		float wetL = spectrumBiquadProcess(&ext->lp, harmL, 0);
		float wetR = spectrumBiquadProcess(&ext->lp, harmR, 1);
		// ViPER topology is parallel-add: wet augments the original, not dry/wet crossfade.
		float outL = (inL * ext->dryMix + wetL * ext->wetMix) * ext->postGain;
		float outR = (inR * ext->dryMix + wetR * ext->wetMix) * ext->postGain;
		if (ext->safetyEnabled)
		{
			outL = tanhf(outL * 1.5f) * (1.0f / 1.5f);
			outR = tanhf(outR * 1.5f) * (1.0f / 1.5f);
		}
		if (!isfinite(outL)) outL = 0.0f;
		if (!isfinite(outR)) outR = 0.0f;
		jdsp->tmpBuffer[0][i] = outL;
		jdsp->tmpBuffer[1][i] = outR;
	}
}
