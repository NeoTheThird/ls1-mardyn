//
// Created by Kruegener on 8/29/2018.
//

#ifndef MARDYN_TRUNK_VELOCITYABSPROFILE_H
#define MARDYN_TRUNK_VELOCITYABSPROFILE_H

#include "ProfileBase.h"
#include "../KartesianProfile.h"

/**
 * @brief Outputs the magnitude of the velocity per bin specified by Sampling grid in KartesianProfile.
 */
class VelocityAbsProfile : public ProfileBase {
public:
    ~VelocityAbsProfile() final = default;
    void record(Molecule& mol, unsigned long uID) final  {
        double absV = 0.0;
        double v;
        for(unsigned short d = 0; d < 3; d++){
            v = mol.v(d);
            absV += v*v;
        }
        absV = sqrt(absV);
        _localProfile[uID] += absV;
    }
    void collectAppend(DomainDecompBase *domainDecomp, unsigned long uID) final {
        domainDecomp->collCommAppendLongDouble(_localProfile[uID]);
    }
    void collectRetrieve(DomainDecompBase *domainDecomp, unsigned long uID) final {
        _globalProfile[uID] = domainDecomp->collCommGetLongDouble();
    }
    void output(string prefix) final;
    void reset(unsigned long uID) final  {
        _localProfile[uID] = 0.0;
        _globalProfile[uID] = 0.0;
    }
    // set correct number of communications needed for this profile
    int comms() final {return 1;}
};

#endif //MARDYN_TRUNK_VELOCITYABSPROFILE_H
