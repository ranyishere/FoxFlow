#ifndef DGGML_Microtubule_TYPES_HPP
#define DGGML_Microtubule_TYPES_HPP
#include "YAGL_Graph.hpp" 
#include "YAGL_Node.hpp" 
#include "SpatialData3D.hpp" 
namespace Microtubule {
struct Type {
	template <class Archive>
	void serialize(Archive& archive) {
	}
};
struct StartType {
	float start_location[3];
	template <class Archive>
	void serialize(Archive& archive) {
		archive(start_location);
	}
};
struct Boundary {
	float boundary_location[3];
	template <class Archive>
	void serialize(Archive& archive) {
		archive(boundary_location);
	}
};
	struct FractureSegment : Type {
		double fflow_6b3cfc[3];
		double fflow_122ce5;
		double fflow_6e2a05[3];
		void* operator[](std::size_t index) const {
			if (index == 0) return (void*)&fflow_6b3cfc;
			if (index == 1) return (void*)&fflow_122ce5;
			if (index == 2) return (void*)&fflow_6e2a05;
			throw std::out_of_range("Index out of bounds");
		};

	template <class Archive>
	void serialize(Archive& archive) {
		archive(fflow_6b3cfc);
		archive(fflow_122ce5);
		archive(fflow_6e2a05);
	}

};	struct FractureSegmentEnd : Type {
		double fflow_8227fa[3];
		double fflow_e957c5;
		double fflow_30aaee[3];
		void* operator[](std::size_t index) const {
			if (index == 0) return (void*)&fflow_8227fa;
			if (index == 1) return (void*)&fflow_e957c5;
			if (index == 2) return (void*)&fflow_30aaee;
			throw std::out_of_range("Index out of bounds");
		};

	template <class Archive>
	void serialize(Archive& archive) {
		archive(fflow_8227fa);
		archive(fflow_e957c5);
		archive(fflow_30aaee);
	}

};	struct Junction : Type {
		double fflow_390a7c[3];
		double fflow_b78c6a;
		double fflow_f00b58[3];
		void* operator[](std::size_t index) const {
			if (index == 0) return (void*)&fflow_390a7c;
			if (index == 1) return (void*)&fflow_b78c6a;
			if (index == 2) return (void*)&fflow_f00b58;
			throw std::out_of_range("Index out of bounds");
		};

	template <class Archive>
	void serialize(Archive& archive) {
		archive(fflow_390a7c);
		archive(fflow_b78c6a);
		archive(fflow_f00b58);
	}

};	struct PressureSegment : Type {
		double fflow_26d4e5[3];
		double fflow_37b7e1;
		void* operator[](std::size_t index) const {
			if (index == 0) return (void*)&fflow_26d4e5;
			if (index == 1) return (void*)&fflow_37b7e1;
			throw std::out_of_range("Index out of bounds");
		};

	template <class Archive>
	void serialize(Archive& archive) {
		archive(fflow_26d4e5);
		archive(fflow_37b7e1);
	}

};	struct PressureSource : Type {
		double fflow_583de4[3];
		double fflow_f25522;
		double fflow_297c72[3];
		void* operator[](std::size_t index) const {
			if (index == 0) return (void*)&fflow_583de4;
			if (index == 1) return (void*)&fflow_f25522;
			if (index == 2) return (void*)&fflow_297c72;
			throw std::out_of_range("Index out of bounds");
		};

	template <class Archive>
	void serialize(Archive& archive) {
		archive(fflow_583de4);
		archive(fflow_f25522);
		archive(fflow_297c72);
	}

};	struct RockStart : Type {
		double fflow_a9c368[3];
		double fflow_a13208;
		double fflow_992eaa;
		double fflow_d20285;
		void* operator[](std::size_t index) const {
			if (index == 0) return (void*)&fflow_a9c368;
			if (index == 1) return (void*)&fflow_a13208;
			if (index == 2) return (void*)&fflow_992eaa;
			if (index == 3) return (void*)&fflow_d20285;
			throw std::out_of_range("Index out of bounds");
		};

	template <class Archive>
	void serialize(Archive& archive) {
		archive(fflow_a9c368);
		archive(fflow_a13208);
		archive(fflow_992eaa);
		archive(fflow_d20285);
	}

};	struct Rock : Type {
		double fflow_e73169[3];
		double fflow_9ca0e1;
		void* operator[](std::size_t index) const {
			if (index == 0) return (void*)&fflow_e73169;
			if (index == 1) return (void*)&fflow_9ca0e1;
			throw std::out_of_range("Index out of bounds");
		};

	template <class Archive>
	void serialize(Archive& archive) {
		archive(fflow_e73169);
		archive(fflow_9ca0e1);
	}

};	using graph_type = YAGL::Graph<std::size_t,	SpatialNode3D<StartType,Boundary,FractureSegment,FractureSegmentEnd,Junction,PressureSegment,PressureSource,RockStart,Rock>>;
};
#endif