/*
 *    Copyright (c) 2015-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_ACTIONS_H_
#define SRC_INCLUDE_ACTIONS_H_

#include <algorithm>
#include <stdexcept>
#include <utility>
#include <vector>

#include "action.h"
#include "forwarddeclarations.h"

namespace smash {

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
  /// Default constructor, creating an empty Actions object.
  Actions() {}
  /**
   * Creates a new Actions object from an ActionList.
   *
   * The actions are stored in a heap and not sorted. The entries of
   * the ActionList are rendered invalid by this constructor.
   *
   * \param[in] action_list The ActionList from which to construct the Actions
   *                    object
   */
  explicit Actions(ActionList&& action_list) : data_(std::move(action_list)) {
    std::make_heap(data_.begin(), data_.end(), cmp);
  }

  /// Cannot be copied
  Actions(const Actions&) = delete;
  /// Cannot be copied
  Actions& operator=(const Actions&) = delete;

  /// \return whether the list of actions is empty.
  bool is_empty() const { return data_.empty(); }

  /**
   * Return the first action in the list and removes it from the list.
   *
   * \throw RuntimeError if the list is empty.
   */
  ActionPtr pop() {
    if (data_.empty()) {
      throw std::runtime_error("Empty actions list!");
    }
    std::pop_heap(data_.begin(), data_.end(), cmp);
    ActionPtr act = std::move(data_.back());
    data_.pop_back();
    return act;
  }

  /**
   * Insert a list of actions into this object.
   *
   * They're inserted at the right places to keep the complete list a heap.
   *
   * \param[in] new_acts The actions that will be inserted.
   */
  void insert(ActionList&& new_acts) {
    for (auto& a : new_acts) {
      insert(std::move(a));
    }
  }

  /**
   * Insert an action into this container.
   *
   * The function makes sure that the action is inserted at the right place.
   *
   * \param[in] action The action to insert.
   */
  void insert(ActionPtr&& action) {
    data_.push_back(std::move(action));
    std::push_heap(data_.begin(), data_.end(), cmp);
  }

  /// \return Number of actions.
  ActionList::size_type size() const { return data_.size(); }

  /// Delete all actions.
  void clear() { data_.clear(); }

  /// \return an iterator to the earliest action.
  std::vector<ActionPtr>::const_reverse_iterator begin() const {
    return data_.crbegin();
  }

  /// \return an iterator to the place following the last action.
  std::vector<ActionPtr>::const_reverse_iterator end() const {
    return data_.crend();
  }

 private:
  /**
   * Compare two action pointer such that the maximum is the most recent
   * action.
   *
   * \param[in] a First action
   * \param[in] b Second action
   * \return Whether the first action will be executed later than the second.
   */
  static bool cmp(const ActionPtr& a, const ActionPtr& b) {
    return a->time_of_execution() > b->time_of_execution();
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

}  // namespace smash

#endif  // SRC_INCLUDE_ACTIONS_H_
