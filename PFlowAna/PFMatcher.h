/*
 * PFMatcher.h
 *
 *  Created on: 03.04.2014
 *      Author: tlodd
 */

#ifndef PFMATCHER_H_
#define PFMATCHER_H_

#include "PFMatchInterfaces.h"
#include "Objects.h"

class DistanceProvider;

namespace PFMatch {

struct MatchDistance {
  MatchDistance(ICluster* cluster, double distance, bool isMatch):
    _cluster(cluster), _distance(distance), _isMatch(isMatch) { }
  ICluster* _cluster;
  double _distance;
  bool _isMatch;
};

class TrackClusterMatcher {
public:
  TrackClusterMatcher(ITrackPositionProvider* trackPosition,
                      IClusterPositionProvider* clusterPosition,
                      IDistanceCalculator* distanceCalculator,
                      double matchCut);
  virtual ~TrackClusterMatcher();

  MatchDistance match(ITrack* track, ICluster* cluster);

  template<class ClusterType>
  MatchDistance bestMatch(ITrack* track, const std::vector<ClusterType*>& clusters);
  template<class ClusterType>
  std::vector<MatchDistance> allMatches(ITrack* track, const std::vector<ClusterType*>& clusters);

private:
  DistanceProvider* _distanceProvider;
  double _matchCut;
};

template<class ClusterType>
MatchDistance TrackClusterMatcher::bestMatch(ITrack* track, const std::vector<ClusterType*>& clusters) {
  ClusterType* bestCluster(0);
  double bestDistance(_matchCut);
  unsigned int nClusters(clusters.size());
  for (unsigned int iCluster = 0; iCluster < nClusters; ++iCluster){
    ClusterType* thisCluster = clusters[iCluster];
    double thisDistance(_distanceProvider->distanceBetween(track, thisCluster));
    std::cout<<"index ="<< iCluster<< " eta="<< clusters[iCluster]->eta() << " phi="<< clusters[iCluster]->phi()<<"  DeltaR="<<thisDistance<<std::endl;
    if (thisDistance < bestDistance) {
      bestDistance = thisDistance;
      bestCluster  = thisCluster;
    }
  }
  return MatchDistance(bestCluster, bestDistance, bestDistance<_matchCut);
}

template<class ClusterType>
std::vector<MatchDistance> TrackClusterMatcher::allMatches(ITrack* track, const std::vector<ClusterType*>& clusters) {
  std::vector<MatchDistance> result;

  unsigned int nClusters(clusters.size());
  for (unsigned int iCluster = 0; iCluster < nClusters; ++iCluster){

    ClusterType* thisCluster = clusters[iCluster];
    double thisDistance(_distanceProvider->distanceBetween(track, thisCluster));

    if (thisDistance < _matchCut) {
      result.push_back(MatchDistance(thisCluster, thisDistance, true));
    }
  }

  return result;
}


//class PlainDRMatcher: public TrackClusterMatcher {
//  PlainDRMatcher(double matchCut) :
//      TrackClusterMatcher(new ClusterPlainEtaPhiProvider(),
//                          new TrackEtaPhiInFixedLayersProvider(EMB2, EME2),
//                          new EtaPhiSqDistanceCalculator(),
//                          matchCut) {
//  }
//  virtual ~PlainDRMatcher() { }
//};
//
// Or simply:
// TrackClusterMatcher* plainDRMatcher = TrackClusterMatcher(new ClusterPlainEtaPhiProvider(),
//                                                           new TrackEtaPhiInFixedLayersProvider(EMB2, EME2),
//                                                           new EtaPhiSqDistanceCalculator(),
//                                                           matchCut);


}

#endif /* PFMATCHER_H_ */
