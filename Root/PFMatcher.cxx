/*
 * PFMatcher.cxx
 *
 *  Created on: 03.04.2014
 *      Author: tlodd
 */

#include <PFlowAna/PFMatcher.h>

namespace PFMatch {
  
  TrackClusterMatcher::TrackClusterMatcher(ITrackPositionProvider* trackPosition,
					   IClusterPositionProvider* clusterPosition,
					   IDistanceCalculator* distanceCalculator, double matchCut) :
    _distanceProvider(new DistanceProvider(trackPosition, clusterPosition, distanceCalculator)),
    _matchCut(matchCut) { }
  
  TrackClusterMatcher::~TrackClusterMatcher() { delete _distanceProvider; }
  
  MatchDistance TrackClusterMatcher::match(ITrack* track, ICluster* cluster) {
    double distance = _distanceProvider->distanceBetween(track, cluster);
    return MatchDistance(cluster, distance, distance<_matchCut);
  }
  
}
