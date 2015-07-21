/*
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_ACTIONS_H_
#define SRC_INCLUDE_ACTIONS_H_

#include <algorithm>
#include <stdexcept>
#include <vector>

#include "action.h"
#include "forwarddeclarations.h"
#include "action.h"

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
  /** Default constructor, creating an empty Actions object. */
  Actions() {}
  /**
   * Creates a new Actions object from an ActionList.
   *
   * The actions are sorted before they are stored in the object. The entries of
   * the ActionList are rendered invalid by this constructor.
   *
   * \param action_list The ActionList from which to construct the Actions
   *                    object
   */
  explicit Actions(ActionList&& action_list)
      : data_(std::move(action_list)) {
    sort(data_);
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
   * Returns the first action in the list and removes it from the list.
   *
   * Throws runtime_error if the list is empty.
   */
  ActionPtr pop() {
    if (data_.empty()) {
      throw std::runtime_error("Empty actions list!");
    }
    ActionPtr act = std::move(data_.back());
    data_.pop_back();
    return std::move(act);
  }

  /**
   * Insert a list of actions into this object.
   *
   * They're inserted at the right places to keep the complete list sorted.
   *
   * \param new_acts The actions that will be inserted.
   */
  void insert(ActionList&& new_acts) {
    if (new_acts.empty()) {
      return;
    }
    sort(new_acts);

    const size_t old_end = data_.size();
    data_.insert(data_.end(), std::make_move_iterator(new_acts.begin()),
                 std::make_move_iterator(new_acts.end()));

    // merge the two lists while preserving the ordering
    std::inplace_merge(data_.begin(), data_.begin() + old_end, data_.end(),
                       [](const ActionPtr &a,
                          const ActionPtr &b) { return *b < *a; });
  }

 private:
  /**
   * Sort the actions such that the first action is at the end of the list.
   *
   * \param action_list The list to sort.
   */
  static void sort(ActionList& action_list) {
    std::sort(action_list.begin(), action_list.end(),
          [](const ActionPtr &a, const ActionPtr &b) { return *b < *a; });
  }

  /**
   * Dynamic data.
   *
   * Vector is likely the best container type here. Because std::sort requires
   * random access iterators. Any linked data structure (e.g. list) thus
   * requires a less efficient sort algorithm.
   */
  std::vector<ActionPtr> data_;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_ACTIONS_H_
