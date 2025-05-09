/*
    nucleus_engine.hpp  -  Don Cross <cosinekitty@gmail.com>

    https://github.com/cosinekitty/sapphire
*/

#pragma once
#include <algorithm>
#include "sapphire_engine.hpp"

namespace Sapphire
{
    inline PhysicsVector EffectiveVelocity(const PhysicsVector& vel, float speedLimit)
    {
        float rawSpeed = Magnitude(vel);
        if (rawSpeed < speedLimit / 1.0e+6)
            return vel;     // prevent division by zero, or by very small numbers

        float effSpeed = BicubicLimiter(rawSpeed, speedLimit);
        return (effSpeed/rawSpeed) * vel;
    }

    struct Particle
    {
        PhysicsVector pos;
        PhysicsVector vel;
        PhysicsVector force;
        float mass = 1.0e-3f;

        bool isFinite() const
        {
            // We don't check `force` because it is a temporary part of the calculation.
            // Any problems in `force` will show up in `pos` and `vel`.
            return pos.isFinite3d() && vel.isFinite3d();
        }
    };


    using NucleusDcRejectFilter = StagedFilter<float, 3>;
    const float DefaultCornerFrequencyHz = 30;


    class NucleusEngine
    {
    private:
        const float max_dt = 0.0005f;
        std::vector<Particle> curr;
        std::vector<Particle> next;
        float magneticCoupling = 0.0f;
        float speedLimit = 1000.0f;
        AutomaticGainLimiter agc;
        bool enableAgc = false;
        int fixedOversample = 0;                // 0 = calculate oversample, >0 = specify oversampling count
        std::vector<float> outputBuffer;        // allows feeding output data through the Automatic Gain Limiter.
        float aetherSpin = 0;
        float aetherVisc = 0;

        // DC reject state (consider moving into a separate class...)
        bool enableDcReject = false;
        const int crossfadeLimit = 8000;
        int crossfadeCounter{};
        float mixFilt{};
        std::vector<NucleusDcRejectFilter> filterArray;     // 3 filters per particle: (x, y, z)
        bool filtersNeedReset = false;

        // Fast inverse square root approximation (Quake III / Carmack technique)
        inline float fastInvSqrt(float number) 
        {
            const float threehalfs = 1.5f;
            float x2 = number * 0.5f;
            float y = number;
            
            // Evil floating point bit level hacking
            int i = *(int*)&y;
            i = 0x5f3759df - (i >> 1); // Magic number for approximation
            y = *(float*)&i;
            
            // Newton-Raphson iteration for accuracy
            y = y * (threehalfs - (x2 * y * y)); // One iteration
            
            return y;
        }
        
        // Fast square root approximation
        inline float fastSqrt(float x) 
        {
            return x * fastInvSqrt(x);
        }

        void calculateForces(std::vector<Particle>& array)
        {
            const float overlapDistance2 = 1.0e-8f;  // Square of overlapDistance for direct comparison
            const int n = static_cast<int>(numParticles());

            // Reset all forces to zero using direct field access for better performance
            for (int i = 0; i < n; ++i)
            {
                array[i].force[0] = 0.0f;
                array[i].force[1] = 0.0f;
                array[i].force[2] = 0.0f;
            }

            // Special case for 0 or 1 particles (no forces to calculate)
            if (n < 2) return;

            // Fast path if no magnetic coupling (purely electrostatic forces)
            if (magneticCoupling == 0.0f)
            {
                for (int i = 0; i+1 < n; ++i)
                {
                    Particle& a = array[i];
                    const float ax = a.pos[0];
                    const float ay = a.pos[1];
                    const float az = a.pos[2];
                    
                    for (int j = i+1; j < n; ++j)
                    {
                        Particle& b = array[j];
                        
                        // Calculate displacement vector components
                        const float dx = b.pos[0] - ax;
                        const float dy = b.pos[1] - ay;
                        const float dz = b.pos[2] - az;
                        
                        // Calculate distance squared directly
                        const float dist2 = dx*dx + dy*dy + dz*dz;

                        // Skip calculation if particles are too close
                        if (dist2 > overlapDistance2)
                        {
                            // Use fast approximation for microcontrollers
                            const float invDist = fastInvSqrt(dist2);
                            const float dist = 1.0f / invDist;
                            const float invDist3 = invDist * invDist * invDist;
                            
                            // Force calculation factor: (dist - 1/dist³)
                            const float electroFactor = dist - invDist3;
                            
                            // Calculate and accumulate forces
                            const float fx = electroFactor * dx;
                            const float fy = electroFactor * dy;
                            const float fz = electroFactor * dz;
                            
                            // Apply forces directly with scalar operations
                            a.force[0] += fx;
                            a.force[1] += fy;
                            a.force[2] += fz;
                            
                            b.force[0] -= fx;
                            b.force[1] -= fy;
                            b.force[2] -= fz;
                        }
                    }
                }
                return;
            }

            // Regular path with magnetic coupling
            for (int i = 0; i+1 < n; ++i)
            {
                Particle& a = array[i];
                const float ax = a.pos[0];
                const float ay = a.pos[1];
                const float az = a.pos[2];
                
                // Precompute and cache effective velocity components
                PhysicsVector av_temp = EffectiveVelocity(a.vel, speedLimit);
                const float avx = av_temp[0];
                const float avy = av_temp[1];
                const float avz = av_temp[2];
                
                for (int j = i+1; j < n; ++j)
                {
                    Particle& b = array[j];
                    
                    // Calculate displacement vector components
                    const float dx = b.pos[0] - ax;
                    const float dy = b.pos[1] - ay;
                    const float dz = b.pos[2] - az;
                    
                    // Calculate distance squared directly
                    const float dist2 = dx*dx + dy*dy + dz*dz;

                    // Skip calculation if particles are too close
                    if (dist2 > overlapDistance2)
                    {
                        // Use fast inverse square root approximation
                        const float invDist = fastInvSqrt(dist2);
                        const float dist = 1.0f / invDist;
                        const float invDist3 = invDist * invDist * invDist;
                        
                        // Force calculation factor: (dist - 1/dist³)
                        const float electroFactor = dist - invDist3;
                        
                        // Cache velocity calculation results (only calculated when needed)
                        PhysicsVector bv_temp = EffectiveVelocity(b.vel, speedLimit);
                        const float dvx = bv_temp[0] - avx;
                        const float dvy = bv_temp[1] - avy;
                        const float dvz = bv_temp[2] - avz;
                        
                        // Calculate magnetic component efficiently
                        const float magFactor = magneticCoupling * invDist3;
                        
                        // Compute force components with minimal operations
                        const float fx = electroFactor * dx + magFactor * (dvy*dz - dvz*dy);
                        const float fy = electroFactor * dy + magFactor * (dvz*dx - dvx*dz);
                        const float fz = electroFactor * dz + magFactor * (dvx*dy - dvy*dx);
                        
                        // Apply forces
                        a.force[0] += fx;
                        a.force[1] += fy;
                        a.force[2] += fz;
                        
                        b.force[0] -= fx;
                        b.force[1] -= fy;
                        b.force[2] -= fz;
                    }
                }
            }
        }

        void extrapolate(float dt)
        {
            const int n = static_cast<int>(numParticles());

            for (int i = 0; i < n; ++i)
            {
                const Particle& p1 = curr.at(i);
                Particle& p2 = next.at(i);

                // F = m*a  ==>  a = F/m.
                PhysicsVector acc = p1.force / p1.mass;

                // Estimate the net velocity change over the interval.
                PhysicsVector dV = dt * acc;

                // Calculate new position using mean velocity change over the interval.
                // We must use the bicubic limiter / "effective velocity" to avoid explosions.

                PhysicsVector v2 = p1.vel + dV/2;
                p2.pos = p1.pos + (dt * EffectiveVelocity(v2, speedLimit));

                // Calculate the velocity at the end of the time interval.
                p2.vel = p1.vel + dV;
            }
        }

        void step(float dt, float friction)
        {
            const int n = static_cast<int>(numParticles());

            // Calculate the forces at the existing configuration.
            calculateForces(curr);

            // Do a naive extrapolation to the midpoint of the time interval.
            // We assume the resulting configuration closely approximates the mean
            // conditions over the whole time interval.
            extrapolate(dt / 2);

            // Calculate forces after moving halfway.
            calculateForces(next);

            // Pretend like the midpoint forces apply at the beginning of the time interval.
            for (int i = 0; i < n; ++i)
                curr[i].force = next[i].force;

            // Extrapolate to the full time interval.
            extrapolate(dt);

            // Update the current state to the calcuted next state.
            // Apply friction at the same time.
            for (int i = 0; i < n; ++i)
            {
                curr[i] = next[i];
                curr[i].vel *= friction;
            }
        }

        float filter(float sampleRate, int i, int k, float x)
        {
            #ifndef METAMODULE
            if (mixFilt > 0)
            {
                // DC rejection is enabled, or we are crossfading.
                NucleusDcRejectFilter& filt = filterArray.at(3*i + k);
                float y = filtersNeedReset ? filt.SnapHiPass(x) : filt.UpdateHiPass(x, sampleRate);
                return (1-mixFilt)*x + mixFilt*y;
            }
            #endif
            return x;
        }

    public:
        explicit NucleusEngine(std::size_t _nParticles)
            : curr(_nParticles)
            , next(_nParticles)
            , outputBuffer(3 * _nParticles)     // (x, y, z) position vectors
            , filterArray(3 * _nParticles)
        {
            initialize();
        }

        void initialize()
        {
            crossfadeCounter = 0;
            enableFixedOversample(1);
            setAgcEnabled(true);
            setDcRejectEnabled(true);
            setDcRejectCornerFrequency(DefaultCornerFrequencyHz);
            setAetherSpin();
            setAetherVisc();
            filtersNeedReset = true;     // anti-click measure: eliminate step function being fed through filters!

            // The caller is responsible for resetting particle states.
            // For example, the caller might want to call SetMinimumEnergy(engine) after calling this function.
        }

        void setDcRejectCornerFrequency(float cutoffFrequencyHz)
        {
            for (NucleusDcRejectFilter& f : filterArray)
                f.SetCutoffFrequency(cutoffFrequencyHz);
        }

        void resetAfterCrash()      // called when infinite/NAN output is detected, to pop back into the finite world
        {
            filtersNeedReset = true;
            agc.initialize();

            const int n = static_cast<int>(numParticles());
            for (int i = 0; i < n; ++i)
                for (int k = 0; k < 3; ++k)
                    output(i, k) = 0;

            // The caller is responsible for resetting particle states.
            // For example, the caller might want to call SetMinimumEnergy(engine) after calling this function.
        }

        void enableFixedOversample(int n)
        {
            fixedOversample = std::max(1, n);
        }

        void enableAutomaticOversample()
        {
            fixedOversample = 0;
        }

        int getOversamplingRate() const
        {
            return fixedOversample;
        }

        bool getAgcEnabled() const
        {
            return enableAgc;
        }

        void setAgcEnabled(bool enable)
        {
            if (enable && !enableAgc)
            {
                // If the AGC isn't enabled, and caller wants to enable it,
                // re-initialize the AGC so it forgets any previous level it had settled on.
                agc.initialize();
            }
            enableAgc = enable;
        }

        void setAgcLevel(float level)
        {
            agc.setCeiling(level);
        }

        float getAetherSpin() const
        {
            return aetherSpin;
        }

        void setAetherSpin(float s = 0)
        {
            aetherSpin = std::clamp(s, -1.0f, +1.0f);
        }

        float getAetherVisc() const
        {
            return aetherVisc;
        }

        void setAetherVisc(float v = 0)
        {
            aetherVisc = std::clamp(v, 0.0f, +1.0f);
        }

        double getAgcDistortion() const     // returns 0 when no distortion, or a positive value correlated with AGC distortion
        {
            return enableAgc ? (agc.getFollower() - 1.0) : 0.0;
        }

        bool getDcRejectEnabled() const
        {
            return enableDcReject;
        }

        void setDcRejectEnabled(bool enable)
        {
            if (enable != enableDcReject)
            {
                // Trigger a cross-fade, to prevent clicking in audio streams.
                enableDcReject = enable;
                crossfadeCounter = crossfadeLimit;

                // Force resetting filters when we process the first input sample.
                // This will make them "snap" to the initial DC state.
                if (enable)
                    filtersNeedReset = true;
            }
        }

        void setMagneticCoupling(float mc)
        {
            magneticCoupling = mc;
        }

        void update(float dt, float halflife, float sampleRate, float gain)
        {
            // Use oversampling to keep the time increment within stability limits.
            // Allow the caller to specify the exact oversampling rate, or we allow
            // the caller to let us figure it out for them (variable performance though).
            int n = fixedOversample;
            if (n < 1)
            {
                // Automatic adjustment of oversampling is enabled.
                n = static_cast<int>(std::ceil(dt / max_dt));
                if (n < 1) n = 1;   // should never happen, but be careful
            }

            // Iterate the model over each oversampled step.
            const double et = dt / n;
            const float friction = std::pow(0.5, static_cast<double>(et)/halflife);
            for (int i = 0; i < n; ++i)
                step(et, friction);

            // Prepare for any crossfading betweeen raw signals and DC rej signals.
            if (enableDcReject || (crossfadeCounter > 0))
            {
                // We will mix each DC-reject signal with the corresponding raw signal using a linear crossfade.
                mixFilt = static_cast<float>(crossfadeCounter) / static_cast<float>(crossfadeLimit);
                if (enableDcReject)
                    mixFilt = 1 - mixFilt;      // toggle the crossfade direction

                if (crossfadeCounter > 0)
                    --crossfadeCounter;
            }

            // Copy outputs, and apply optional DC rejection.
            // For a couple hundred samples after toggling the DC rejection option,
            // there is a linear crossfade between the raw signal and the filtered
            // signal, to reduce audio popping.
            const int nparticles = static_cast<int>(curr.size());
            for (int i = 0; i < nparticles; ++i)
            {
                const Particle& p = curr.at(i);
                for (int k = 0; k < 3; ++k)
                    output(i, k) = gain * filter(sampleRate, i, k, p.pos[k]);
            }

            if (mixFilt > 0)
                filtersNeedReset = false;       // the `filter` function has finished resetting all the filters by now

            if (enableAgc)
                agc.process(sampleRate, outputBuffer);
        }

        float& output(int p, int k)
        {
            return outputBuffer.at(3*p + k);
        }

        std::size_t numParticles() const
        {
            return curr.size();
        }

        Particle& particle(int index)
        {
            return curr.at(index);
        }

        const Particle& particle(int index) const
        {
            return curr.at(index);
        }
    };
}
