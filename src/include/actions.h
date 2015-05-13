/*
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_ACTIONS_H_
#define SRC_INCLUDE_ACTIONS_H_

#include <algorithm>
#include <forward_list>

#include "forwarddeclarations.h"

namespace Smash {

/**
 * \ingroup data
 *
 * The Actions class abstracts the storage and manipulation of actions.
 *
 * \note
 * The Actions object cannot be copied, because it does not make sense
 * semantically. Move semantics make sense and can be implemented when needed.
 */
class Actions {
 public:
  /**
   * Creates a new Actions object from an ActionList.
   *
   * The actions are sorted before they are stored in the object. The entries of
   * the ActionList are rendered invalid by this constructor.
   *
   * \param action_list The ActionList from which to construct the Actions
   *                    object
   */
  Actions(ActionList&& action_list, float current_time) {
    // sort the actions
    std::sort(action_list.begin(), action_list.end(),
          [](const ActionPtr &a, const ActionPtr &b) { return *a < *b; });

    // move the ActionPtrs from action_list to data_
    // and make the internal time of the actions global
    for (auto it = action_list.rbegin(); it != action_list.rend(); ++it) {
      (*it)->make_time_global(current_time);
      data_.push_front(std::move(*it));
    }
  }

  /**
   * Cannot be copied
   */
  Actions(const Actions&) = delete;
  /**
   * Cannot be copied
   */
  Actions& operator=(const Actions&) = delete;

  /**
   * Returns whether the list of actions is empty.
   */
  bool is_empty() const { return data_.empty(); }

  /**
   * Returns the first action in the list.
   */
  ActionPtr pop() {
    auto it = data_.begin();
    if (it == data_.end()) {
      return nullptr;
    }
    ActionPtr act = std::move(*it);
    data_.pop_front();
    return std::move(act);
  }

  /**
   * Insert a list of actions into this object.
   *
   * They're inserted at the right places to keep the complete list sorted.
   *
   * \param new_acts The actions that will be inserted.
   */
  void insert(ActionList&& new_acts, float current_time) {
    if (new_acts.empty()) {
      // nothing to do
      return;
    }
    // sort first
    std::sort(new_acts.begin(), new_acts.end(),
          [](const ActionPtr &a, const ActionPtr &b) { return *a < *b; });

    // correctly set the time of execution
    for (const auto &act : new_acts) {
      act->make_time_global(current_time);
    }

    // iterator that points at before the beginning
    // this is necessary because there is only a function insert_after
    // for which we need the element before the place where we want to insert
    auto before_it0 = data_.before_begin();

    // first action of the list
    auto new_it = new_acts.begin();
    float new_time = (*new_it)->time_of_execution();

    for (auto it0 = data_.begin(); it0 != data_.end(); ++it0) {
      if (new_time < (*it0)->time_of_execution()) {
        // the action at it0 is after our action -> insert it before it0
        it0 = data_.insert_after(before_it0, std::move(*new_it));
        new_acts.erase(new_it);
        // check if there are actions left
        if (new_acts.empty()) {
          return;
        }
        // look at the next action of the list
        new_it = new_acts.begin();
        new_time = (*new_it)->time_of_execution();
      }
      ++before_it0;
    }
    // insert the rest
    for (new_it = new_acts.begin(); new_it != new_acts.end(); ++new_it) {
      data_.insert_after(before_it0, std::move(*new_it));
      ++before_it0;
    }
  }

 private:
  /**
   * Dynamic data.
   */
  std::forward_list<ActionPtr> data_;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_ACTIONS_H_
